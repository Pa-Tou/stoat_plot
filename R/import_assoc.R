#' Import association results from STOAT
#'
#' @description Read the TSV (potentially gzipped) file produced by STOAT and load it.
#' @param assoc_file Path to STOAT's output (*assoc.pvalues.tsv.gz)
#'
#' @return A polished data.frame with a subset of the TSV file.
#' @name import_assoc
#' @export

import_assoc <- function(assoc_file) {

  ## don't load all columns (add more if/when we need them)
  ## all_cols = c('#CHR', 'START_OFFSET', 'END_OFFSET',
  ##              'START_NODE', 'END_NODE',
  ##              'P_FISHER', 'P_CHI2', 'ALLELE_COUNT_PER_PHENO',
  ##              'ALLELE_LENGTHS', 'GENE', 'P',
  ##              'ALLELE_COUNT', 'DEPTH')
  sel_cols = c('#CHR'='character',
               'START_OFFSET'='numeric', 
               'START_NODE'='character',
               'END_NODE'='character',
               'P_FISHER'='numeric',
               'P_CHI2'='numeric',
               'ALLELE_COUNT_PER_PHENO'='character',
               'ALLELE_LENGTHS'='character',
               'GENE'='character',
               'P'='numeric',
               'ALLELE_COUNT'='character')

  ## read a few rows to check the header
  df = readr::read_delim(assoc_file, n_max=3, progress=FALSE, show_col_types=FALSE)
  sel_cols = sel_cols[which(names(sel_cols) %in% colnames(df))]
  
  ## read the full file
  df = readr::read_delim(assoc_file, col_select=names(sel_cols),
                         ## n_max =100000,
                         col_types=sel_cols, show_col_types=FALSE)

  ## Remove leading # in headers (if present)
  colnames(df)[which(colnames(df) == '#CHR')] = 'CHR'

  return(df)
}
