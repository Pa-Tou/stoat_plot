#' P-value Histogram for STOAT GWAS Results
#' @description Generate histogram of P-values from STOAT GWAS results TSV.
#'
#' @param input Path to input TSV file.
#' @param output Filename to save the output plot (default: "pvalue_distribution_plot.png").
#' @param p_threshold Maximum P-value threshold to include in the plot (default: 0.1).
#' @param bin Number of bins in the histogram (default: 200).
#' @param column_names Column name to use for p-values (default: ""). If empty, will use "P" or "P_CHI2" if available.
#'
#' @return Saves a histogram plot as an image file.
#' @name plot_pvalue_hist
#' @export
plot_pvalue_hist <- function(input, column_names = "", p_threshold = 0.1, bin = 200, output = "pvalue_distribution_plot.png") {

  # Read data
  data <- read.table(input, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                     check.names = FALSE, comment.char = "")
  colnames(data)[1] <- sub("^#", "", colnames(data)[1])  # remove leading # from first column name

  # Determine the column to use for p-values
  if (column_names != "") {
    if (!(column_names %in% colnames(data))) {
      stop(paste("Column", column_names, "not found in the input file."))
    }
    p_col <- column_names
  } else {
    if ("P" %in% colnames(data)) {
      p_col <- "P"
    } else if ("P_CHI2" %in% colnames(data)) {
      p_col <- "P_CHI2"
    } else {
      stop("Neither 'P' nor 'P_CHI2' column found in the input file.")
    }
  }

  # Convert p-value column to numeric
  data[[p_col]] <- as.numeric(data[[p_col]])

  if (any(is.na(data[[p_col]]))) {
    warning("NA values detected in the p-value column. They will be excluded from the plot.")
  }

  # Filter data by removing NA and applying threshold
  data_filtered <- data[!is.na(data[[p_col]]) & data[[p_col]] <= p_threshold, ]

  if (nrow(data_filtered) == 0) {
    stop("No valid P-values found within the specified threshold.")
  }

  # Plot histogram using ggplot2
  p <- ggplot2::ggplot(data_filtered, ggplot2::aes_string(x = p_col)) +
    ggplot2::geom_histogram(bins = bin, fill = "skyblue", color = "black") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Distribution of P-values (0 -", p_threshold, ")"),
      x = "P-value",
      y = "Frequency"
    )

  # Save plot
  ggplot2::ggsave(output, plot = p, width = 8, height = 6, dpi = 300)
}
