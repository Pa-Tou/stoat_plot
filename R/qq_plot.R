#' Q-Q Plot for GWAS Results
#' @description Generate QQ plot from STOAT GWAS results using P or P_CHI2 column.
#'
#' @param input Path to the input TSV file (must contain a column named 'P').
#' @param output Filename for the output PNG plot (default: "qq_plot.png").
#' @param column_names Column name to use for p-values (default: ""). If empty, will use "P" or "P_CHI2" if available.
#'
#' @return Saves a Q-Q plot image.
#' @name qq_plot
#' @export
qq_plot <- function(input, column_names="", output = "qq_plot.png") {

  # Read input file
  data <- read.table(input, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "" )
  colnames(data)[1] <- sub("^#", "", colnames(data)[1])  # remove leading #

  # Determine the column to use for p-values
  if (column_names != "") {
    if (!(column_names %in% colnames(data))) {
      stop(paste("Column", column_names, "not found in the input file."))
    }
    p_column <- column_names
  } else {
    if ("P" %in% colnames(data)) {
      p_column <- "P"
    } else if ("P_CHI2" %in% colnames(data)) {
      p_column <- "P_CHI2"
    } else {
      stop("Neither 'P' nor 'P_CHI2' column found in the input file.")
    }
  }

  # Convert p-value column to numeric
  pvals <- as.numeric(data[[p_column]])

  if (any(is.na(pvals) & pvals < 0 & pvals > 1)) {
    warning("Invalide values detected in the p-value column. They will be excluded from the plot.")
  }

  # Clean P-values
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
  pvals <- pmax(pvals, 1e-300)  # Avoid log(0)

  # Expected and observed -log10(P)
  n <- length(pvals)
  expected <- -log10((1:n) / (n + 1))
  observed <- -log10(sort(pvals))

  # Data frame for plotting
  plot_df <- data.frame(Expected = expected, Observed = observed)

  # Plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Expected, y = Observed)) +
    ggplot2::geom_point(size = 1.5, color = "steelblue") +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    ggplot2::labs(
      title = "QQ Plot",
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  # Save plot
  ggplot2::ggsave(output, plot = p, width = 6, height = 6)
}
