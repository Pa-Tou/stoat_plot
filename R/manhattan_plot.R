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
  ## Input sanity checks
  ## ---------------------------
  if (!is.null(start) && !is.null(end) && start > end) {
    stop("start must be <= end")
  }

  # return error if gwas_file does not exist
  if (!file.exists(gwas_file)) {
    stop("gwas_file does not exist: ", gwas_file)
  }

  # return an error if p_threshold is not a positive number
  if (!is.numeric(p_threshold) || p_threshold <= 0) {
    stop("p_threshold must be a positive number")
  }

  # -----------------------------
  # Read header of GWAS file
  # and validate format
  # -----------------------------
  con <- file(gwas_file, open = "r")
  on.exit(close(con))

  header_line <- NULL
  chr_list <- character(0)
  line_count <- 0

  repeat {
    line <- readLines(con, n = 1)

    if (length(line) == 0) {
      stop("Unexpected end of file before finding '#START_NODE'.")
    }

    line_count <- line_count + 1

    # Capture column header
    if (grepl("^#START_NODE\\b", line)) {
      header_line <- line
      break
    }

    # Capture chromosome section
    if (grepl("^#REFS\\b", line)) {

      repeat {
        line <- readLines(con, n = 1)

        if (length(line) == 0) {
          stop("Unexpected end of file inside '#REFS' section.")
        }

        line_count <- line_count + 1

        # End of chromosome section
        if (grepl("^#SNARLS\\b", line)) {
          break
        }

        # Remove leading '#' and store chromosome name
        chr_name <- sub("^#", "", line)
        chr_list[length(chr_list) + 1] <- chr_name
      }
    }
  }

  if (is.null(header_line)) {
    stop("Column header '#START_NODE' not found in GWAS file.")
  }

  # Remove leading '#'
  header_line <- sub("^#", "", header_line)

  # Split columns
  col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]

  # Expected fixed header
  expected_cols <- c(
    "START_NODE",
    "END_NODE",
    "REF_INDEX",
    "START_OFFSET",
    "END_OFFSET",
    "DEPTH",
    "ALLELE_LENGTHS",
    "WALKS",
    "SEQUENCES",
    "P"
  )

  # Check exact match
  if (!identical(col_names[1:9], expected_cols)) {
    stop(
      paste0(
        "Invalid gwas file format.\n",
        "Expected first 10 columns:\n",
        paste(expected_cols, collapse = "\t")
      )
    )
  }

  gwas_data <- read.table(
    gwas_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names,
    skip = line_count,
    comment.char = "",
    quote = "",
    fill = FALSE,
    check.names = FALSE
  )

  # -----------------------------
  # Filter by chromosome if provided
  # -----------------------------
  if (!is.null(chr)) {

    # Find index(es) of requested chromosome(s)
    chr_idx <- match(chr, chr_list)

    if (any(is.na(chr_idx))) {
      stop(sprintf(
        "Chromosome(s) not found in header: %s",
        paste(chr[is.na(chr_idx)], collapse = ", ")
      ))
    }

    # Convert to 0-based index (since REF_INDEX is 0-based)
    chr_idx <- chr_idx - 1

    # Filter data
    gwas_data <- gwas_data[gwas_data$REF_INDEX %in% chr_idx, , drop = FALSE]

    if (nrow(gwas_data) == 0) {
      stop("No variants found for the specified chromosome(s).")
    }
  } else {
    # If no chromosome specified, add CHR column based on REF_INDEX
    gwas_data$CHR <- factor(
      chr_list[gwas_data$REF_INDEX + 1],
      levels = chr_list
    )
  }

  # -----------------------------
  # Prepare gwas_data
  # -----------------------------
  gwas_data$START_OFFSET <- as.integer(gwas_data$START_OFFSET)
  gwas_data$P <- pmax(as.numeric(gwas_data[[p_column]]), 1e-300)

  gwas_data <- gwas_data[!is.na(gwas_data$START_OFFSET) & !is.na(gwas_data$P), ]

  if (!is.null(start)) {
    gwas_data <- gwas_data[gwas_data$START_OFFSET >= start, ]
  }

  if (!is.null(end)) {
    gwas_data <- gwas_data[gwas_data$START_OFFSET <= end, ]
  }

  # -----------------------------
  # Prepare plotting gwas_data
  # -----------------------------
  gwas_data <- data.frame(
    CHR = gwas_data$CHR,
    BP = gwas_data$START_OFFSET,
    P = gwas_data$P,
    stringsAsFactors = FALSE
  )

  gwas_data$logp <- -log10(gwas_data$P)

  # -----------------------------
  # X-axis handling
  # -----------------------------
  if (!is.null(chr)) {
    gwas_data$xpos <- gwas_data$BP
    axis_df <- NULL
    x_label <- paste0(chr, " position (bp)")
  } else {
    gwas_data$CHR <- factor(gwas_data$CHR, levels = unique(gwas_data$CHR))
    gwas_data <- gwas_data[order(gwas_data$CHR, gwas_data$BP), ]

    chr_lengths <- tapply(gwas_data$BP, gwas_data$CHR, max)
    chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
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

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(gwas_data, aes(x = gwas_data$xpos, y = gwas_data$logp)) +
    geom_point(
      aes(color = gwas_data$CHR),
      alpha = 0.6,
      size = 0.7
    ) +
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
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)
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
