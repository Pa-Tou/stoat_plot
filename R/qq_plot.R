#' Q-Q Plot for GWAS Results
#' @description Generate QQ plot from STOAT GWAS results using P or P_CHI2 column.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline
#' @importFrom utils read.table
#'
#' @param gwas_file Path to the output stoat GWAS TSV file.
#' @param p_column Column name to use for p-values (default: ""). If empty, will use "P" or "P_CHI2" if available.
#' @param output Filename for the output PNG plot (default: "qq_plot.png").
#'
#' @return Saves a Q-Q plot image.
#' @name qq_plot
#' @export

qq_plot <- function(gwas_file, 
              p_column = "P", 
              output = "qq_plot.png") {

  ## ---------------------------
  ## Input checks
  ## ---------------------------
  if (!file.exists(gwas_file)) {
    stop("gwas_file does not exist: ", gwas_file)
  }

  ## ---------------------------
  ## Read GWAS file
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

  # Check p-value column
  if (!(p_column %in% colnames(gwas_data))) {
    stop(paste("Column", p_column, "not found in the gwas_data file."))
  }

  # Convert p-values
  pvals <- suppressWarnings(as.numeric(gwas_data[[p_column]]))

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
  chisq <- qchisq(1 - pvals, df = 1, lower.tail = FALSE)
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
