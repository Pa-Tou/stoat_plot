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

manhattan_plot <- function(input, column_names="", p_threshold=1e-5, output = "manhattan_plot.png") {

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

  pvals <- as.numeric(data[[p_column]])

  # Required columns
  if (!all(c("CHR", "START_POS") %in% colnames(data))) {
    stop("Input file must contain columns: 'CHR' and 'START_POS'")
  }

  # Clean START_POS
  data$START_POS <- as.integer(gsub(",", "", as.character(data$START_POS)))

  # Subset and clean
  data <- data[!is.na(data$CHR) & !is.na(data$START_POS) & !is.na(data[[p_column]]), ]
  data <- data.frame(
    CHR = data$CHR,
    BP = data$START_POS,
    P = pmax(pvals, 1e-300),
    stringsAsFactors = FALSE
  )

  # Ensure chromosome is treated as a factor to control plotting order
  data$CHR <- factor(data$CHR, levels = unique(data$CHR[order(as.character(data$CHR))]))

  # Sort and calculate cumulative START_POSition
  data <- data[order(data$CHR, data$BP), ]
  chr_lengths <- tapply(data$BP, data$CHR, max)
  chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
  names(chr_offsets) <- names(chr_lengths)
  data$cum_bp <- data$BP + chr_offsets[as.character(data$CHR)]

  # Calculate -log10(P)
  data$logp <- -log10(data$P)

  # Midpoints for axis labels
  axis_df <- aggregate(cum_bp ~ CHR, data = data, FUN = function(x) (min(x) + max(x)) / 2)

  logp_threshold <- -log10(p_threshold)

  # Plot
  p <- ggplot(data, aes(x = cum_bp, y = logp)) +
    geom_point(aes(color = CHR), alpha = 0.6, size = 0.7) +
    scale_color_manual(values = rep(c("black", "grey50"), length.out = length(levels(data$CHR)))) +
    scale_x_continuous(breaks = axis_df$cum_bp, labels = axis_df$CHR) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](P)),
      title = "Manhattan Plot"
    ) +
    # Add threshold line
    geom_hline(yintercept = logp_threshold, color = "red", linetype = "dashed") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",  # fixed typo: was "legend.START_POSition"
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)
    )

  # Save plot
  ggsave(output, plot = p, width = 12, height = 4)
}
