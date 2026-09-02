#' Summary of STOAT's results
#'
#' @description Print a few summary metrics from the association results
#' @importFrom stats p.adjust median qchisq setNames
#' @importFrom rlang .data
#'
#' @param assoc Either the data.frame imported by *import_assoc*, or the path to STOAT's output (*assoc.pvalues.tsv.gz)
#'
#' @return In addition to printing the summary metrics, returns a list with the values
#' 
#' @name summary_stoat
#' @examples
#' # prepare the filename
#' assoc_file = system.file('extdata/stoat.quantitative.assoc.pvalues.tsv.gz', package='StoatPlot')
#'
#' # print a summary
#' sum.l = summary_stoat(assoc_file)
#' sum.l 
#' @export

summary_stoat <- function(assoc) {

  ## potentially load the association results
  if(is.character(assoc) && length(assoc) == 1) {
    assoc = import_assoc(assoc)
  }

  ## which column to use for the pvalue?
  pv_col = ifelse(any(colnames(assoc) == 'P'), 'P', 'P_CHI2')
  
  ## correct for multiple test correction
  if(all(colnames(assoc) != 'P_BH')){
    assoc$P_BH = p.adjust(assoc[,pv_col, TRUE], method = "BH")
  }

  ## Genomic inflation factor
  chisq_stat <- median(qchisq(1 - assoc[,pv_col, TRUE], df = 1), na.rm = TRUE)
  inf.lambda <- chisq_stat / qchisq(0.5, df = 1)

  ## compute some statistics on the pvalues
  sum.df = tibble::tibble(Metric="Total variants", Value=nrow(assoc))
  if(any(is.na(assoc[,pv_col]))) {
    sum.df = rbind(
      sum.df,
      tibble::tibble(Metric="Total variants with non-NA PV", Value=sum(!is.na(assoc[,pv_col]))))
  }
  sum.df = rbind(
    sum.df,
    tibble::tibble(Metric="Variant PV<0.01", Value=sum(assoc[,pv_col] < .01, na.rm=TRUE)),
    tibble::tibble(Metric="Variant adjusted PV<0.01", Value=sum(assoc$P_BH < .01, na.rm=TRUE))
  )
  if(any(colnames(assoc) == 'GENE')) {
    sum.df = rbind(
      sum.df,
      tibble::tibble(Metric="Genes PV<0.01", Value=length(unique(assoc$GENE[which(assoc[,pv_col] < .01)]))),
      tibble::tibble(Metric="Genes adjusted PV<0.01", Value=length(unique(assoc$GENE[which(assoc$P_BH < .01)])))
    )
  }
  sum.df = rbind(
    sum.df,
    tibble::tibble(Metric="Genomic inflation factor", Value=inf.lambda)
  )

  ## print as Markdown tables
  cat('\n\n')
  cat(knitr::kable(sum.df, format.args=list(digits=3, scientific=FALSE,
                                            big.mark=',', drop0trailing=TRUE)), sep='\n')

  ## assign a variant type based on PATH_LENGTHS
  ## TODO: move to import_assoc?
  sum.per.type = NULL
  if(any(colnames(assoc) == 'ALLELE_LENGTHS')) {
    assoc$Type = assign_type_from_lengths(assoc$ALLELE_LENGTHS)

    sum.per.type = dplyr::summarize(dplyr::group_by(assoc, .data$Type),
                                    Total=dplyr::n(),
                                    Pv_BH_below_0.01=sum(.data$P_BH < .01))
    sum.per.type$Type = factor(sum.per.type$Type,
                               levels=c("SNP", "MNP", "SV"))
    sum.per.type = sum.per.type[order(sum.per.type$Type),]
    
    ## print as Markdown tables
    cat('\n\n')
    cat(knitr::kable(sum.per.type, format.args=list(big.mark=',')), sep='\n')
    cat('\n\n')
  }
  
  return(list(all=sum.df, per.type=sum.per.type))
}
