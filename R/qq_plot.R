#' Q-Q Plot for GWAS Results
#' @description Generate QQ plot from STOAT GWAS results using P or P_CHI2 column.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline
#' @importFrom utils read.table
#'
#' @param input Path to the input TSV file (must contain a column named 'P').
#' @param column_names Column name to use for p-values (default: ""). If empty, will use "P" or "P_CHI2" if available.
#' @param output Filename for the output PNG plot (default: "qq_plot.png").
#'
#' @return Saves a Q-Q plot image.
#' @name qq_plot
#' @export

qq_plot <- function(input, column_names = "P", output = "qq_plot.png") {

  # Read all lines
  lines <- readLines(input)

  # Identify the header line (starts with #START_NODE)
  header_idx <- grep("^#START_NODE", lines)
  if (length(header_idx) == 0) {
    stop("Header line '#START_NODE' not found in the file.")
  }

  # Extract header and clean '#'
  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]

  # Extract data lines (after header, not starting with #)
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  # Read data into data.frame
  data <- read.table(
    text = data_lines,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names
  )

  # Check p-value column
  if (!(column_names %in% colnames(data))) {
    stop(paste("Column", column_names, "not found in the input file."))
  }

  # Convert p-values
  pvals <- suppressWarnings(as.numeric(data[[column_names]]))

  # Filter valid p-values
  valid <- !is.na(pvals) & pvals > 0 & pvals <= 1
  if (any(!valid)) {
    warning("Invalid p-values detected and removed.")
  }
  pvals <- pvals[valid]
  pvals <- pmax(pvals, 1e-300)

  # Expected and observed -log10(P)
  n <- length(pvals)
  expected <- -log10((1:n) / (n + 1))
  observed <- -log10(sort(pvals))

  plot_df <- data.frame(Expected = expected, Observed = observed)

  # Genomic inflation factor
  chisq <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

  # Plot
  p <- ggplot(plot_df, aes(plot_df$Expected, plot_df$Observed)) +
    geom_point(
      size = 1.5,
      color = "cadetblue3",
      alpha = 0.7
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "red"
    ) +
    labs(
      title = "QQ Plot",
      subtitle = paste0("Genomic inflation factor = ", round(lambda, 3)),
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank()
    )

  ggsave(output, plot = p, width = 6, height = 6)
}
