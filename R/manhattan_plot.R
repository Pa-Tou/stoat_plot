#' Manhattan Plot for GWAS Results
#' @description Generate Manhattan plots from STOAT GWAS results (keeps CHR names like 'chr1', 'chrX', etc.)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_color_manual scale_x_continuous labs theme_bw theme element_blank element_text ggsave
#' @importFrom utils read.table head write.table
#' @importFrom stats aggregate
#'
#' @param gwas_file Path to the output stoat GWAS TSV file.
#' @param p_column Column name to use for p-values (default: "P").
#' @param chr Optional column name for chromosome (default: NULL, will try "CHR").
#' @param start Optional column name for start positions (default: NULL, will try "START_OFFSET").
#' @param end Optional column name for end positions (default: NULL, will try "END_OFFSET").
#' @param p_threshold P-value threshold for the horizontal significance line (default: 1e-5).
#' @param output Path to save the output plot image.
#'
#' @return Saves a Manhattan plot to the specified file.
#' @name manhattan_plot
#' @export

manhattan_plot <- function(gwas_file,
                           p_column = "P",
                           chr = NULL,
                           start = NULL,
                           end = NULL,
                           p_threshold = 1e-5,
                           output = "manhattan_plot.png") {

  ## ---------------------------
  ## Input checks
  ## ---------------------------
  if (!file.exists(gwas_file)) {
    stop("gwas_file does not exist: ", gwas_file)
  }

  if (!is.numeric(p_threshold) || p_threshold <= 0) {
    stop("p_threshold must be a positive number")
  }

  if (!is.null(start) && !is.null(end) && start > end) {
    stop("start must be <= end")
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

  ## ---------------------------
  ## Validate required columns
  ## ---------------------------
  required_cols <- c("CHR", "START_OFFSET", p_column)

  missing_cols <- setdiff(required_cols, colnames(gwas_data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  ## ---------------------------
  ## Prepare columns
  ## ---------------------------
  gwas_data$BP <- as.integer(gwas_data$START_OFFSET)
  gwas_data$P  <- pmax(as.numeric(gwas_data[[p_column]]), 1e-300)
  gwas_data <- gwas_data[!is.na(gwas_data$BP) & !is.na(gwas_data$P), ]

  ## ---------------------------
  ## Filters
  ## ---------------------------
  if (!is.null(chr)) {
    gwas_data <- gwas_data[gwas_data$CHR %in% chr, ]

    if (nrow(gwas_data) == 0) {
      stop("No variants found for the specified chromosome(s).")
    }
  }

  if (!is.null(start)) {
    gwas_data <- gwas_data[gwas_data$BP >= start, ]
  }

  if (!is.null(end)) {
    gwas_data <- gwas_data[gwas_data$BP <= end, ]
  }

  ## ---------------------------
  ## Prepare plotting data
  ## ---------------------------
  gwas_data$logp <- -log10(gwas_data$P)

  ## ---------------------------
  ## X-axis handling
  ## ---------------------------
  if (!is.null(chr)) {

    gwas_data$xpos <- gwas_data$BP
    axis_df <- NULL
    x_label <- paste0(chr, " position (bp)")

  } else {

    gwas_data$CHR <- factor(gwas_data$CHR, levels = unique(gwas_data$CHR))
    gwas_data <- gwas_data[order(gwas_data$CHR, gwas_data$BP), ]

    chr_lengths <- tapply(gwas_data$BP, gwas_data$CHR, max)
    chr_offsets <- c(0, cumsum(chr_lengths)[-length(chr_lengths)])
    names(chr_offsets) <- names(chr_lengths)

    gwas_data$xpos <- gwas_data$BP + chr_offsets[as.character(gwas_data$CHR)]

    axis_df <- aggregate(
      xpos ~ CHR,
      data = gwas_data,
      FUN = function(x) mean(range(x))
    )

    x_label <- "Chromosome"
  }

  logp_threshold <- -log10(p_threshold)

  ## ---------------------------
  ## Plot
  ## ---------------------------
  p <- ggplot(gwas_data, aes(x = xpos, y = logp)) +
    geom_point(aes(color = CHR), alpha = 0.6, size = 0.7) +
    geom_hline(
      yintercept = logp_threshold,
      color = "red",
      linetype = "dashed"
    ) +
    labs(
      x = x_label,
      y = expression(-log[10](P)),
      title = "Manhattan Plot"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  if (is.null(chr)) {
    p <- p +
      scale_x_continuous(
        breaks = axis_df$xpos,
        labels = axis_df$CHR
      ) +
      scale_color_manual(
        values = rep(
          c("cadetblue3", "darkcyan"),
          length.out = length(levels(gwas_data$CHR))
        )
      )
  }

  ggsave(output, plot = p, width = 12, height = 4)
}
