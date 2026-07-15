#' Q-Q plot of the p-values
#' @description Generate QQ plot of the observed vs expected pvalues
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline
#' @importFrom rlang .data
#'
#' @param assoc Either the data.frame imported by *import_assoc*, or the path to STOAT's output (*assoc.pvalues.tsv.gz)
#' @param output_file If not NULL, the name of the output image where to save the plot (image type guessed from the file name).
#'
#' @return a ggplot object. Saves a file too if output_file is provided.
#' @name qq_plot
#' @examples
#' # prepare the filename
#' assoc_file = system.file('extdata/stoat.quantitative.assoc.pvalues.tsv.gz', package='StoatPlot')
#'
#' # make the QQ-plot
#' qq_plot(assoc_file)
#' @export

qq_plot <- function(assoc, output_file=NULL) {

  ## potentially load the association results
  if(is.character(assoc) && length(assoc) == 1) {
    assoc = import_assoc(assoc)
  }

  ## which column to use for the pvalue?
  pv_col = ifelse(any(colnames(assoc) == 'P'), 'P', 'P_CHI2')

  ## Genomic inflation factor
  chisq_stat <- median(qchisq(1 - assoc[, pv_col, TRUE], df = 1), na.rm = TRUE)
  lambda <- chisq_stat / qchisq(0.5, df = 1)
    
  ## handle null pvalues
  null_pvs = which(assoc[, pv_col] == 0)
  if (length(null_pvs) > 0) {
    min.pv = min(assoc[-null_pvs, pv_col])
    min.pv = 10^(floor(log10(min.pv) - 5))
    warning(length(null_pvs), ' p-values changed from 0 to ', min.pv)
    assoc[null_pvs, pv_col] = min.pv
  }
  
  # expected and observed -log10(P)
  df <- tibble::tibble(expected=stats::ppoints(nrow(assoc)),
                       observed=sort(assoc[, pv_col, TRUE]))

  ## if very large, downsample pvalues above .01
  total_variants = nrow(df)
  if(total_variants > 1e6) {
    n_downsampled = round(total_variants/5e4)
    df = df[which(df$expected < .0001 | 1:nrow(df) %% n_downsampled == 0), ]
  }
  
  ## plot
  ggp <- ggplot(df, aes(x=-log10(.data$expected), y=-log10(.data$observed))) +
    geom_point(size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(x = "expected -log10(P)", y = "observed -log10(P)",
         caption = paste0("Genomic inflation factor = ", round(lambda, 3))) +
    theme_bw()

  if(!is.null(output_file)) {
    ggsave(output_file, plot=ggp, width=6, height=6)
  }

  return(ggp)
}
