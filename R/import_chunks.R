##' Import a slice of the snarl genotype file (sorted, bgzipped and
##' indexed with Tabix). For each snarl, a table with sample and
##' allele counts is return in a list, along with other information
##' about the snarl.
##'
##' The list returned as one element per snarl. Each snarl element has
##' a *GT* data.frame with a *sample* column and a column for each
##' alelle (in the form al0, al1, etc). Each row in this *GT*
##' data.frame correspond to a sample and inform each allele's
##' count. That table could be combined with the phenotype and used
##' for a regression test, for example.
##' @title Import a chunk of the snarl genotypes for a queried region
##' @param genotype_file path to the snarl genotype file (sorted and
##'   indexed with tabix)
##' @param chrom chromosome name
##' @param start_offset start position
##' @param end_offset end position
##' @param keep_nonvariant should we keep variants with only one
##'   allele present? Default FALSE.
##' @param long_format should the table we returned in long form. Default FALSE.
##' @return a list of snarls with some information (boundary nodes)
##'   and the genotype table
##' @export
import_genotype_chunk <- function(genotype_file, chrom, start_offset, end_offset,
                                  keep_nonvariant=FALSE, long_format=FALSE) {
  ## input sanity checks
  if (!file.exists(genotype_file)) {
    stop("genotype_file does not exist: ", genotype_file)
  }

  ## set up tabix file connection
  gt.tbx <- Rsamtools::TabixFile(genotype_file)

  ## read the header to get the sequence IDs
  first.id = NULL
  last.id = NULL
  gt.h = NULL
  read.n = 100
  while(length(last.id) == 0) {
    gt.h = scan(genotype_file, '', sep='\n', n=read.n, quiet=TRUE)
    first.id = grep('#REFS', gt.h) + 1
    last.id = grep('#SNARLS', gt.h)
    read.n = read.n * 2
  }
  seqs.ids = 1:(last.id-first.id) - 1
  names(seqs.ids) = gsub('^#', '', gt.h[first.id:(last.id-1)])

  ## header of the snarl genotype part
  gt.snarls.h = unlist(strsplit(gt.h[last.id+1], '\t'))
  
  ## get a slice using tabix
  query.gr = GenomicRanges::GRanges(seqs.ids[chrom],
                                    IRanges::IRanges(start_offset, end_offset))
  res <- strsplit(Rsamtools::scanTabix(gt.tbx, param=query.gr)[[1]], '\t')

  ## prepare sample vector from the header
  samples = gsub('(.+)#.+', '\\1', gt.snarls.h[10:length(gt.snarls.h)])
  
  ## format the result into a list with some info and the genotype data.frame
  res.ll = lapply(res, function(rr) {
    ## start a long-format table with the allele presence
    alleles = rr[10:length(gt.snarls.h)]
    df = tibble::tibble(sample=samples,
                        allele=as.numeric(ifelse(alleles == '.', NA, alleles)))

    ## skip non-variants
    if(!keep_nonvariant && length(unique(df$allele)) == 1) {
      return(NULL)
    }
    
    ## aggregate into allele counts
    df = dplyr::summarize(dplyr::group_by(df, .data$sample, .data$allele),
                          count=dplyr::n(), .groups='drop')

    ## turn to wide-format table
    if (!long_format) {
      df = tidyr::pivot_wider(df, names_from=.data$allele, names_prefix='al',
                              values_from=.data$count, values_fill=0)
    }
    
    ## remove the NA/missing allele count
    if(any(colnames(df) == 'alNA')) {
      df$alNA = NULL
    }
    
    list(START_NODE=rr[1], END_NODE=rr[2], WALKS=rr[grep('WALKS', gt.snarls.h)],
         ALLELE_LENGTHS=rr[grep('ALLELE_LENGTHS', gt.snarls.h)],
         GT=df)
  })

  ## don't return NULL elements (snarls with one allele)
  res.null = sapply(res.ll, is.null)
  if(length(res.ll) == 0 || all(res.null)) {
    # if all snarls are NULL or the list was empty
    return(list())
  }
  if(any(res.null)) {
    res.ll = res.ll[which(!res.null)]
  }
  
  return(res.ll)
}

##' Combine a phenotype and potentially some covariables with the
##' allele count information from import_genotype_chunk. If paths are
##' given to phenotypr or covariables, will try to read the files.
##' @title Combine the allele counts for snarl with phenotype and covariables
##' @param snarl_ac data.frame with the allele counts for the snarl
##' @param phenotype Either a data.frame with the phenotype for each sample, or the path to the phenotype file
##' @param covariables Either a data.frame with the covariables for each sample, or the path to the covariables file
##' @param gene_name gene of interest, if the phenotype is a matrix of gene expression.
##' @return a data.frame with the new phenotype and covariables columns
##' @export
combine_chunk_phenotype_covars <- function(snarl_ac, phenotype, covariables=NULL,
                                           gene_name=NULL) {
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
    if (is.null(gene_name)) {
      stop("Phenotype file error: must specify a gene name with gene_name.")
    }
    if (all(gene_name != phenotype$gene_name)) {
      stop("Phenotype file error: no gene name found (", gene_name, ")")
    }
    
    query.pheno = phenotype[which(phenotype$gene_name == gene_name),]
    query.pheno$gene_name = NULL
    query.pheno = tibble::tibble(sample=colnames(query.pheno),
                                 phenotype=as.numeric(query.pheno))
  }

  ## merge phenotype
  df = merge(snarl_ac, query.pheno)

  ## eventually merge covariables
  if (!is.null(covariables)){
    ## potentially load the covariables information
    if(is.character(covariables) && length(covariables) == 1) {
      if (!file.exists(covariables)) {
        stop("covariable file does not exist: ", covariables)
      }
      covariables = readr::read_tsv(covariables, progress=FALSE,
                                    show_col_types=FALSE)
    }

    ## add a prefix to the covariables, for convenience?
    colnames(covariables) = ifelse(colnames(covariables) == 'SAMPLE',
                                   'sample', paste0('covar', colnames(covariables)))
    
    ## combine with current table
    df = merge(df, covariables)
  }

  return(df)
}


##' Prepare the phenotype/genotype/covariate table before a test
##'
##' Removes the sample name column, adds a column for the total allele
##' count of this snarl, removes duplicated columns, and removes one
##' allele. The total allele count is to avoid the effect in the
##' parent snarl to leak to the current snarl.
##' @param test_table data.frame generated by combine_chunk_phenotype_covars
##' @return a data.frame ready for a test
##' @export
preprocess_test_table <- function(test_table) {
  ## all variable names except the phenotype and sample IDs
  var.names = setdiff(colnames(test_table), c('sample', 'phenotype'))

  ## find constant columns
  col.var = sapply(var.names, function(varn) {
    stats::var(as.numeric(test_table[[varn]]))
  })
  var.to.rm = var.names[col.var == 0]
  var.names = var.names[col.var != 0]

  ## add total allele count (if it varies)
  al.names = grep('^al', var.names, value=TRUE)
  if (length(al.names) == 0){
    ## stop if we've already removed all alleles (unlikely but possible)
    return(data.frame())
  }
  total.ac = rowSums(as.matrix(test_table[,al.names]))
  if (stats::var(total.ac) > 0) {
    test_table$total_allele_count = total.ac
    var.names = c(var.names, 'total_allele_count')
  }
  
  ## find duplicated columns
  col.cat = sapply(var.names, function(varn) {
    paste(as.numeric(test_table[[varn]]), collapse='_')
  })
  var.to.rm = c(var.to.rm, var.names[which(duplicated(col.cat))])
  var.names = var.names[which(!duplicated(col.cat))]

  ## find the first allele
  al.names = grep('^al', var.names, value=TRUE)
  if (length(al.names) == 0){
    ## stop if we've already removed all alleles (unlikely but possible)
    return(data.frame())
  } else if (length(al.names) > 1) {
    ## only remove an allele if there are two?
    var.to.rm = c(var.to.rm, head(al.names, 1))
  }
  
  ## remove the variables
  for(varn in var.to.rm) {
    test_table[[varn]] = NULL
  }

  return(test_table)
}
