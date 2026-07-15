#' P-value Histogram for STOAT GWAS Results
#' @description Generate histogram of P-values from STOAT GWAS results TSV.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline aes_string geom_histogram
#' @importFrom rlang .data
#' @importFrom utils read.table
#'
#' @param gwas_file Path to the output stoat GWAS TSV file.
#' @param output Filename to save the output plot (default: "pvalue_distribution_plot.png").
#' @param min Minimun P-value threshold to include in the plot (default: 0).
#' @param max Maximum P-value threshold to include in the plot (default: 1.0).
#' @param bins Number of bins in the histogram (default: 200).
#' @param p_column Column name to use for p-values (default: ""). If empty, will use "P" or "P_CHI2" if available.
#'
#' @return Saves a histogram plot as an image file.
#' @name plot_pvalue_hist
#' @export

plot_pvalue_hist <- function(gwas_file,
                             p_column = "P",
                             min = 0,
                             max = 1.0,
                             bins = 100,
                             output = NULL) {

  ## ---------------------------
  ## Input checks
  ## ---------------------------
  if (!file.exists(gwas_file)) {
    stop("gwas_file does not exist: ", gwas_file)
  }

  if (!is.numeric(min) || !is.numeric(max) || min < 0 || max > 1 || min >= max) {
    stop("Invalid p-value range. Must satisfy: 0 <= min < max <= 1.")
  }

  ## ---------------------------
  ## Read file
  ## ---------------------------
  gwas_data <- read.table(
    gwas_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )

  colnames(gwas_data) <- sub("^#", "", colnames(gwas_data))

  if (!(p_column %in% colnames(gwas_data))) {
    stop("Column not found: ", p_column)
  }

  ## ---------------------------
  ## Convert p-values
  ## ---------------------------
  gwas_data[[p_column]] <- suppressWarnings(
    as.numeric(gwas_data[[p_column]])
  )

  if (any(is.na(gwas_data[[p_column]]))) {
    warning("NA or non-numeric values detected and removed.")
  }

  ## ---------------------------
  ## Filter valid p-values
  ## ---------------------------
  data_filtered <- gwas_data[
    !is.na(gwas_data[[p_column]]) &
      gwas_data[[p_column]] >= min &
      gwas_data[[p_column]] <= max,
  ]

  if (nrow(data_filtered) == 0) {
    stop("No valid P-values found within the specified range.")
  }

  ## ---------------------------
  ## Plot
  ## ---------------------------
  p <- ggplot2::ggplot(data_filtered,
                       ggplot2::aes(x = .data[[p_column]])) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = "cadetblue3",
      color = "darkcyan",
      alpha = 0.7
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste0("Distribution of P-values (", min, " - ", max, ")"),
      x = "P-value",
      y = "Frequency"
    )

  if(!is.null(output)) {
    ggplot2::ggsave(output, plot = p, width = 8, height = 6, dpi = 300)
  }

  return(p)
}
