#' Manhattan Plot for GWAS Results
#' @description Generate Manhattan plots from STOAT GWAS results (keeps CHR names like 'chr1', 'chrX', etc.)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_color_manual scale_x_continuous labs theme_bw theme element_blank element_text ggsave geom_abline
#' @importFrom utils read.table
#'
#' @param input Path to the input GWAS TSV file.
#' @param output Path to save the output plot image.
#' @param column_names Column name to use for p-values (default: ""). If empty, will use "P" or "P_CHI2" if available.
#' @param p_threshold P-value threshold for the horizontal significance line.
#'
#' @return Saves a Manhattan plot to the specified file.
#' @name manhattan_plot
#' @export

manhattan_plot <- function(input,
                            p_column = "P",
                            chr = NULL,
                            start = NULL,
                            end = NULL,
                            p_threshold = 1e-5,
                            output = "manhattan_plot.png") {

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

  # Read data
  data <- read.table(
    text = data_lines,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names
  )

  # -----------------------------
  # Validate required columns
  # -----------------------------
  if (!(p_column %in% colnames(data))) {
    stop(paste("Column:", p_column, "not found in the input file."))
  }

  if (!("START_OFFSET" %in% colnames(data))) {
    stop("Input file must contain column: 'START_OFFSET'")
  }

  # -----------------------------
  # Prepare data
  # -----------------------------
  data$START_OFFSET <- as.integer(data$START_OFFSET)
  data$P <- pmax(as.numeric(data[[p_column]]), 1e-300)

  data <- data[!is.na(data$START_OFFSET) & !is.na(data$P), ]

  # Assign chromosome
  if (!is.null(chr)) {
    data$CHR <- chr
  } else if ("REF_INDEX" %in% colnames(data)) {
    data$CHR <- paste0("ref", data$REF_INDEX)
  } else {
    stop("Chromosome information missing: provide 'chr' or 'REF_INDEX' column.")
  }

  # -----------------------------
  # Genomic filters
  # -----------------------------
  if (!is.null(start)) {
    data <- data[data$START_OFFSET >= start, ]
  }

  if (!is.null(end)) {
    data <- data[data$START_OFFSET <= end, ]
  }

  if (!is.null(start) && !is.null(end) && start > end) {
    stop("start must be <= end")
  }

  # -----------------------------
  # Prepare plotting data
  # -----------------------------
  data <- data.frame(
    CHR = data$CHR,
    BP = data$START_OFFSET,
    P = data$P,
    stringsAsFactors = FALSE
  )

  data$logp <- -log10(data$P)

  # -----------------------------
  # X-axis handling
  # -----------------------------
  if (!is.null(chr)) {
    data$xpos <- data$BP
    axis_df <- NULL
    x_label <- paste0(chr, " position (bp)")
  } else {
    data$CHR <- factor(data$CHR, levels = unique(data$CHR))
    data <- data[order(data$CHR, data$BP), ]

    chr_lengths <- tapply(data$BP, data$CHR, max)
    chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
    names(chr_offsets) <- names(chr_lengths)

    data$xpos <- data$BP + chr_offsets[as.character(data$CHR)]

    axis_df <- aggregate(
      xpos ~ CHR,
      data = data,
      FUN = function(x) mean(range(x))
    )

    x_label <- "Chromosome"
  }

  logp_threshold <- -log10(p_threshold)

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(data, aes(x = xpos, y = logp)) +
    geom_point(aes(color = CHR), alpha = 0.6, size = 0.7) +
    geom_hline(yintercept = logp_threshold, color = "red", linetype = "dashed") +
    labs(
      x = x_label,
      y = expression(-log[10](P)),
      title = "Manhattan Plot"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)
    )

  if (is.null(chr)) {
    p <- p +
      scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$CHR) +
      scale_color_manual(
        values = rep(c("black", "grey50"),
                     length.out = length(levels(data$CHR)))
      )
  }

  ggsave(output, plot = p, width = 12, height = 4)
}
