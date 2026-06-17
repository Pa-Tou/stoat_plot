#' Manhattan Plot for GWAS Results
#' @description Generate Manhattan plots from STOAT GWAS results (keeps CHR names like 'chr1', 'chrX', etc.)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_text labs theme_minimal ggsave
#' @importFrom stats aggregate
#'
#' @param input Path to the input GWAS TSV file.
#' @param p_column Column name to use for p-values (default: "P").
#' @param chr If provided, plot single chromosome positions instead of concatenating chromosomes.
#' @param start Optional minimum BP to include.
#' @param end Optional maximum BP to include.
#' @param p_threshold P-value threshold for the horizontal significance line (default: 1e-5).
#' @param output Path to save the output plot image.
#' @param interactive If TRUE, return an interactive plotly object (requires plotly in Suggests).
#' @param apply_gc If TRUE, apply genomic control (divide chi-square by lambda) before plotting.
#' @param show_fdr If TRUE, add a horizontal line for BH FDR threshold (alpha = 0.05 by default).
#'
#' @return Invisibly returns ggplot or plotly object and saves static PNG when output provided.
#' @name manhattan_plot
#' @export

manhattan_plot <- function(input,
                           p_column = "P",
                           chr = NULL,
                           start = NULL,
                           end = NULL,
                           p_threshold = 1e-05,
                           output = "manhattan_plot.png",
                           interactive = FALSE,
                           apply_gc = FALSE,
                           show_fdr = TRUE) {
  # Lightweight header parse (header line starts with #CHR)
  lines <- readLines(input)
  header_idx <- grep("^#CHR", lines)
  if (length(header_idx) == 0) stop("Header line '#CHR' not found in the input file.")
  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  # Use data.table::fread for performance on large files
  dt <- data.table::fread(paste(data_lines, collapse = "\n"), sep = "\t", header = FALSE, col.names = col_names, data.table = FALSE)

  if (!(p_column %in% colnames(dt))) stop(paste("Column:", p_column, "not found in the input file."))
  if (!("START_OFFSET" %in% colnames(dt))) stop("Input file must contain column: 'START_OFFSET'")

  dt$START_OFFSET <- as.integer(dt$START_OFFSET)
  dt$P_raw <- pmax(as.numeric(dt[[p_column]]), 1e-300)

  # Compute lambda for genomic control if requested
  chisq <- stats::qchisq(1 - dt$P_raw, df = 1)
  lambda <- median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
  if (apply_gc && !is.na(lambda) && lambda > 0) {
    chisq_gc <- chisq / lambda
    dt$P <- pmax(1 - stats::pchisq(chisq_gc, df = 1), 1e-300)
  } else {
    dt$P <- dt$P_raw
  }

  df <- data.frame(CHR = dt$CHR, BP = dt$START_OFFSET, P = dt$P, stringsAsFactors = FALSE)
  df <- df[!is.na(df$BP) & !is.na(df$P), ]
  df$logp <- -log10(df$P)

  # Handle concatenated x-axis when multiple chromosomes
  if (!is.null(chr)) {
    df$xpos <- df$BP
    x_label <- paste0(chr, " position (bp)")
  } else {
    # use gtools::mixedsort to get human order if available
    if (requireNamespace("gtools", quietly = TRUE)) {
      chr_levels <- gtools::mixedsort(unique(df$CHR))
    } else {
      chr_levels <- sort(unique(df$CHR))
    }
    df$CHR <- factor(df$CHR, levels = chr_levels)
    df <- df[order(df$CHR, df$BP), ]
    chr_lengths <- tapply(df$BP, df$CHR, max)
    chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
    names(chr_offsets) <- names(chr_lengths)
    df$xpos <- df$BP + chr_offsets[as.character(df$CHR)]
    axis_df <- aggregate(xpos ~ CHR, data = df, FUN = function(x) mean(range(x)))
    x_label <- "Position (Mbp)"
  }

  logp_threshold <- -log10(p_threshold)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = xpos / 1e6, y = logp, color = CHR)) +
    ggplot2::geom_point(alpha = 0.6, size = 0.7) +
    ggplot2::geom_hline(yintercept = logp_threshold, color = "red", linetype = "dashed") +
    ggplot2::labs(x = x_label, y = expression(-log[10](P)), title = "Manhattan Plot") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   axis.title = ggplot2::element_text(face = "bold"),
                   axis.text = ggplot2::element_text(color = "black"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none")

  # Add FDR line if requested
  if (show_fdr && nrow(df) > 0) {
    bh <- stats::p.adjust(df$P, method = "BH")
    # find P value corresponding to BH threshold 0.05
    idx <- which(bh <= 0.05)
    if (length(idx) > 0) {
      pval_fdr <- max(df$P[idx], na.rm = TRUE)
      p <- p + ggplot2::geom_hline(yintercept = -log10(pval_fdr), color = "blue", linetype = "dotted")
    }
  }

  # Save static image
  if (!interactive && !is.null(output)) {
    ggplot2::ggsave(output, plot = p, width = 14, height = 6)
    message("Saved plot: ", output)
    invisible(p)
  } else if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) stop("plotly is required for interactive plots; install it or set interactive=FALSE")
    pltly <- plotly::ggplotly(p)
    return(pltly)
  } else {
    invisible(p)
  }
}
