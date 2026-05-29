#!/usr/bin/Rscript
# Run as:
# ./manhattan_plot.R [input file from stoat]  
# Writes file manhattan_plot.png to current directory
# TODO: Add other options as command line options


#' Manhattan Plot for GWAS Results
#' @description Generate Manhattan plots from STOAT GWAS results (keeps CHR names like 'chr1', 'chrX', etc.)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_color_manual scale_x_continuous labs theme_bw theme element_blank element_text ggsave
#' @importFrom utils read.table head write.table
#' @importFrom stats aggregate
#'
#' @param input Path to the input GWAS TSV file.
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

manhattan_plot <- function(input,
                            p_column = "P",
                            chr = NULL,
                            start = NULL,
                            end = NULL,
                            p_threshold = 1.3072e-05,
                            output = "manhattan_plot.png") {

  # -----------------------------
  # Read file lines
  # -----------------------------
  lines <- readLines(input)

  # Detect header line
  header_idx <- grep("^#CHR", lines)
  if (length(header_idx) == 0) {
    stop("Header line '#CHR' not found in the input file.")
  }

  # Parse header
  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]

  # Extract data lines
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  # Read data
  data <- read.delim(
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
  #if (!is.null(chr)) {
  #  data$CHR <- chr
  #} else if ("REF_INDEX" %in% colnames(data)) {
  #  data$CHR <- paste0("ref", data$REF_INDEX)
  #} else {
  #  stop("Chromosome information missing: provide 'chr' or 'REF_INDEX' column.")
  #}

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
    data$CHR <- factor(data$CHR, levels = sort(unique(data$CHR)))
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

    x_label <- "Position (Mbp)"
  }

  logp_threshold <- -log10(p_threshold)

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(data, aes(x = data$BP/1e6, y = data$logp)) +
    geom_point(
      aes(color = data$CHR),
      alpha = 0.6,
      size = 0.7
    ) +
    geom_hline(
      yintercept = logp_threshold,
      color = "red",
      linetype = "dashed"
    ) +
    facet_grid(.~CHR, scales="free", space="free") + 
    labs(
      x = x_label,
      y = expression(-log[10](P)),
      title = "Manhattan Plot"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)
    )

  if (is.null(chr)) {	    p <- p +
      scale_color_manual(
        values = rep(
          c("cadetblue3", "darkcyan"),
          length.out = length(levels(data$CHR))
        )

      )
  }

  ggsave(output, plot = p, width = 20, height = 10)
}

# If this is called as a script, do this. I think this will prevent it from being called when just importing the file
if (sys.nframe() == 0){
    library(tidyverse)
    manhattan_plot(commandArgs(TRUE)[1], output=commandArgs(TRUE)[2])
}
