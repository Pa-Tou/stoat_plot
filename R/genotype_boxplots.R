#' Genotype vs phenotype boxplots
#' @description Generates boxplots to illustrate an association, showing the phenotype distribution for each genotype group.
#'
#' The queried association can be provided either with a list or a
#' data.frame. If a list is provided as input, it must contain the
#' following variables: 'CHR', 'START_OFFSET', 'START_NODE',
#' 'END_NODE'. This information can be found in the TSV files produced
#' by Stoat (either the snarl information or the association
#' results). If the input is a data.frame, only the first row will be
#' used to extract those same information. In practice, this input
#' data.frame could be one row of the association data.frame loaded by
#' *import_assoc*.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_bw ggsave xlab position_jitter geom_label
#' @importFrom rlang .data
#'
#' @param genotype_file path to the snarl genotype file (sorted and indexed with tabix)
#' @param phenotype Either a data.frame with the phenotype for each sample, or the path to the phenotype file
#' @param assoc_query a list or a data frame with one row. See details.
#' @param output_file If not NULL, the name of the output image where to save the plot (image type guessed from the file name).
#' @param by_allele group the samples by allele (an heterozygous sample would be shown twice). Default: FALSE (i.e. grouped by genotype)
#' @param show_n should the sample size be shown at the bottom? Default is TRUE.
#'
#' @return a ggplot object. Saves a file if output_file is provided.
#' @name genotype_boxplots
#' @export

genotype_boxplots <- function(genotype_file, phenotype, assoc_query, output_file=NULL, by_allele=FALSE, show_n=TRUE) {
  
  ## input sanity checks
  if (!file.exists(genotype_file)) {
    stop("genotype_file does not exist: ", genotype_file)
  }

  if (is.list(assoc_query) & !is.data.frame(assoc_query)) {
    ## input is a list, double-check all information is present
    missing.vars = setdiff(c('CHR', 'START_OFFSET', 'START_NODE', 'END_NODE'),
                           names(assoc_query))
    if(length(missing.vars) > 0) {
      stop('Input list (assoc_query) is missing ', paste(missing.vars, collapse=' '))
    }
  }
  if (is.data.frame(assoc_query)) {
    if(nrow(assoc_query) > 1) {
      assoc_query = assoc_query[1,]
      warning('Multiple input associations specified. Using the first one.')
    }
  }

  ## potentially load the phenotype information
  if(is.character(phenotype) && length(phenotype) == 1) {
    if (!file.exists(phenotype)) {
      stop("phenotype_file does not exist: ", phenotype)
    }
    phenotype = readr::read_tsv(phenotype, progress=FALSE, show_col_types=FALSE)
  }

  ## prepare the phenotype if needed
  query.pheno = NULL
  if (any(colnames(phenotype) == 'PHENO')) {
    ## binary or quantitative phenotype
    query.pheno = tibble::tibble(sample=phenotype$SAMPLE, phenotype=phenotype$PHENO)
  } else {
    ## gene expression
    if (all(colnames(phenotype) != 'gene_name')) {
      stop("Phenotype file error: no column gene_name found.")
    }
    query.pheno = phenotype[which(phenotype$gene_name == assoc_query$GENE),]
    query.pheno$gene_name = NULL
    query.pheno = tibble::tibble(sample=colnames(query.pheno),
                                 phenotype=as.numeric(query.pheno))
  }
  

  chunk.l = import_genotype_chunk(genotype_file, assoc_query$CHR,
                                  assoc_query$START_OFFSET, assoc_query$START_OFFSET,
                                  keep_nonvariant=TRUE, long_format=TRUE)
  
  ## find the snarl of interest
  sel.bool = sapply(chunk.l, function(rr) {
    return(rr$START_NODE == assoc_query$START_NODE & rr$END_NODE == assoc_query$END_NODE)
  })
  query.res = chunk.l[[which(sel.bool)]]
  
  ## start a table with the allele for each sample-haplotype
  df = query.res$GT

  ggp = ggplot()
  ## prepare title
  if (any(colnames(assoc_query) == 'GENE')) {
    ggp.title = paste0(assoc_query$GENE, ' - ', assoc_query$CHR, ':', assoc_query$START_OFFSET)
  } else {
    ggp.title = paste0(assoc_query$CHR, ':', assoc_query$START_OFFSET)
  }
  ## prepare allele paths for the caption
  alleles.path = unlist(strsplit(query.res$WALKS, ','))
  alleles.lens = unlist(strsplit(query.res$ALLELE_LENGTHS, ','))
  alleles.path = paste(1:length(alleles.path) - 1, ': ', alleles.path,
                       ' (', alleles.lens, ' bp)', sep='', collapse='\n')
  
  if (by_allele) {
    ## one panel per allele, grouped by allele count

    ## turn to wide-format table
    df = tidyr::pivot_wider(df, names_from=.data$allele,
                            values_from=.data$count, values_fill=0)
    df = tidyr::pivot_longer(df, cols=-.data$sample,
                             names_to='allele', values_to='allele_count')

    ## merge genotypes and phenotype
    df = merge(df, query.pheno)

    ## make graph
    ggp = ggplot(df, aes(x=.data$allele_count, y=.data$phenotype, group=.data$allele_count)) + 
      geom_point(alpha = 0.3, position=position_jitter(.2)) +
      geom_violin(fill = "cadetblue3", alpha = 0.3) +
      geom_boxplot(width = 0.2, outliers=FALSE, outlier.colour = "red",
                   alpha = 0.5, fill = "darkcyan") +
      theme_bw() + labs(caption=alleles.path, title=ggp.title) +
      scale_x_continuous(breaks=0:10) +
      xlab('allele count') + 
      facet_grid(.~allele)

    if(show_n) {
      ## add sample size
      df.n = dplyr::summarize(dplyr::group_by(df, .data$allele, .data$allele_count),
                              allele_count_n=dplyr::n(),
                              .groups='drop')
      df.n$phenotype = min(df$phenotype) - .05 * diff(range(df$phenotype))
      
      ggp = ggp + geom_label(aes(label=.data$allele_count_n), data=df.n, vjust=1, size=3)
    }
    
  } else {
    ## group by genotype
    ## aggregate genotypes
    df = dplyr::summarize(dplyr::group_by(df, .data$sample),
                          allele=paste(sort(rep(.data$allele, .data$count)), collapse='/'), .groups='drop')

    ## merge genotypes and phenotype
    df = merge(df, query.pheno)

    ## make graph
    ggp = ggplot(df, aes(x=.data$allele, y=.data$phenotype, group=.data$allele)) + 
      geom_point(alpha = 0.3, position=position_jitter(.05)) +
      geom_violin(fill = "cadetblue3", alpha = 0.3) +
      geom_boxplot(width = 0.2, outliers=FALSE, outlier.colour = "red",
                   alpha = 0.5, fill = "darkcyan") +
      theme_bw() + labs(caption=alleles.path, title=ggp.title)

    if(show_n) {
      ## add sample size
      df.n = dplyr::summarize(dplyr::group_by(df, .data$allele),
                              allele_n=dplyr::n(),
                              .groups='drop')
      df.n$phenotype = min(df$phenotype) - .05 * diff(range(df$phenotype))
      
      ggp = ggp + geom_label(aes(label=.data$allele_n), data=df.n, vjust=1, size=3)
    }

  }
  
  if(!is.null(output_file)) {
    ggsave(filename=output_file, plot=ggp, width=8, height=6, dpi=300)
  }

  return(ggp)
}
