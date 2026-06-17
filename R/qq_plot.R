#' Q-Q Plot for GWAS Results
#' @description Generate QQ plot from STOAT GWAS results using a p-value column.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal geom_abline ggsave
#'
#' @param input Path to the input TSV file.
#' @param p_column Column name to use for p-values (default: "P_CHI2"). Falls back to "P" if missing.
#' @param output Filename for the output PNG plot (default: "qq_plot.png").
#' @param interactive If TRUE, returns an interactive plotly object (requires plotly in Suggests).
#' @param apply_gc If TRUE, apply genomic control correction to p-values before plotting.
#'
#' @return Saves or returns a Q-Q plot.
#' @name qq_plot
#' @export

qq_plot <- function(input, p_column = "P_CHI2", output = "qq_plot.png", interactive = FALSE, apply_gc = FALSE) {
  lines <- readLines(input)
  header_idx <- grep("^#CHR", lines)
  if (length(header_idx) == 0) stop("Header line '#CHR' not found in the file.")
  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  df <- data.table::fread(paste(data_lines, collapse = "\n"), sep = "\t", header = FALSE, col.names = col_names, data.table = FALSE)

  # choose column
  if (!(p_column %in% colnames(df))) {
    if ("P" %in% colnames(df)) p_column <- "P" else stop(paste("Column", p_column, "not found in the input file."))
  }

  pvals <- suppressWarnings(as.numeric(df[[p_column]]))
  valid <- !is.na(pvals) & pvals > 0 & pvals <= 1
  if (any(!valid)) warning("Invalid p-values detected and removed.")
  pvals <- pvals[valid]
  pvals <- pmax(pvals, 1e-300)

  n <- length(pvals)
  expected <- -log10((1:n) / (n + 1))
  observed <- -log10(sort(pvals))

  # genomic inflation
  chisq <- stats::qchisq(1 - pvals, df = 1)
  lambda <- stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1)

  if (apply_gc && !is.na(lambda) && lambda > 0) {
    chisq_gc <- chisq / lambda
    pvals_gc <- pmax(1 - stats::pchisq(chisq_gc, df = 1), 1e-300)
    observed <- -log10(sort(pvals_gc))
  }

  plot_df <- data.frame(Expected = expected, Observed = observed)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Expected, y = Observed)) +
    ggplot2::geom_point(size = 1.5, color = "cadetblue3", alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "QQ Plot", subtitle = paste0("Genomic inflation factor = ", round(lambda, 3)), x = "Expected -log10(P)", y = "Observed -log10(P)") +
    ggplot2::theme_minimal(base_size = 14)

  if (!interactive) {
    ggplot2::ggsave(output, plot = p, width = 6, height = 6)
    message("Saved plot: ", output)
    invisible(p)
  } else {
    if (!requireNamespace("plotly", quietly = TRUE)) stop("plotly is required for interactive plots; install it or set interactive=FALSE")
    return(plotly::ggplotly(p))
  }
}
