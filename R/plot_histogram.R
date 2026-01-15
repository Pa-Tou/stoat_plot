#' P-value Histogram for STOAT GWAS Results
#' @description Generate histogram of P-values from STOAT GWAS results TSV.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline aes_string geom_histogram
#' @importFrom utils read.table
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

plot_pvalue_hist <- function(input,
                             column_names = "P",
                             p_threshold = 1e-5,
                             bin = 200,
                             output = "pvalue_distribution_plot.png") {

  # -----------------------------
  # Read file lines
  # -----------------------------
  lines <- readLines(input)

  # Detect header line
  header_idx <- grep("^#START_NODE", lines)
  if (length(header_idx) == 0) {
    stop("Header line '#START_NODE' not found in the input file.")
  }

  # Parse header
  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]

  # Extract data lines
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  # Read data table
  data <- read.table(
    text = data_lines,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names
  )

  # -----------------------------
  # Determine p-value column
  # -----------------------------
  if (!(column_names %in% colnames(data))) {
    stop(paste("Column", column_names, "not found in the input file."))
  }
  p_col <- column_names

  # Convert to numeric
  data[[p_col]] <- suppressWarnings(as.numeric(data[[p_col]]))

  if (any(is.na(data[[p_col]]))) {
    warning("NA or non-numeric values detected in the p-value column. They were excluded.")
  }

  # Filter valid p-values
  data_filtered <- data[
    !is.na(data[[p_col]]) &
      data[[p_col]] > 0 &
      data[[p_col]] <= p_threshold,
  ]

  if (nrow(data_filtered) == 0) {
    stop("No valid P-values found within the specified threshold.")
  }

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(data_filtered, aes(x = .data[[p_col]])) +
    geom_histogram(bins = bin, fill = "skyblue", color = "black") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Distribution of P-values (0 –", p_threshold, ")"),
      x = "P-value",
      y = "Frequency"
    )

  # Save plot
  ggsave(output, plot = p, width = 8, height = 6, dpi = 300)
}
