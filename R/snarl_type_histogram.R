
#' Dot Plot of snarl type histogram
#' @description Create a histogram plot from a TSV file containing a `path_length` column.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline geom_bar scale_fill_manual
#' @importFrom utils read.table
#' 
#' @param input Path to the input TSV file.
#' @param output Path to save the output plot image.
#'
#' @return Saves a histogram plot to the specified file.
#' @name snarl_type_histogram
#' @export

snarl_type_histogram <- function(input, output = "snarl_type_histogram.png") {

  con <- file(input, "r")
  on.exit(close(con))

  # --- Find header line ---
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      stop("Unexpected end of file before finding '#START_NODE'.")
    }

    if (grepl("^#START_NODE\\b", line)) {
      header_line <- sub("^#", "", line)
      break
    }
  }

  # Split column names
  col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]

  # Read data from current position
  df <- read.table(
    con,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names,
    comment.char = "",
    quote = "",
    fill = FALSE,
    check.names = FALSE
  )

  if (!"ALLELE_LENGTHS" %in% colnames(df)) {
    stop("Column 'ALLELE_LENGTHS' not found in data.")
  }

  # Classify variant types
  df$Variant_Type <- sapply(df$ALLELE_LENGTHS, function(path) {

    if (is.na(path) || path == ".") {
      return(NA)
    }

    values <- as.numeric(unlist(strsplit(path, "[,/]+")))
    values <- values[!is.na(values)]

    if (length(values) == 0) {
      return(NA)
    }

    max_val <- max(values)

    if (all(values == 1)) {
      return("SNP")
    } else if (max_val <= 50) {
      return("INDEL")
    } else {
      return("SV")
    }
  })

  # Aggregate counts by Variant_Type
  variant_counts <- as.data.frame(table(df$Variant_Type))
  colnames(variant_counts) <- c("Variant_Type", "Count")

  # Ensure all types are present
  all_types <- c("SNP", "INDEL", "SV")
  for (t in all_types) {
    if (!(t %in% variant_counts$Variant_Type)) {
      variant_counts <- rbind(variant_counts, data.frame(Variant_Type = t, Count = 0))
    }
  }

  # Order factor levels
  variant_counts$Variant_Type <- factor(variant_counts$Variant_Type, levels = all_types)

  print(variant_counts)

  # ----------------- PLOT -----------------
plot <- ggplot(
  variant_counts,
  aes(
    x = Variant_Type,
    y = Count,
    fill = Variant_Type
  )
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    alpha = 0.8
  ) +
  labs(
    title = "Snarl Type Distribution",
    x = "Snarl Type",
    y = "Count",
    fill = "Variant Type :"   # <- this fixes legend title
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  scale_fill_manual(
    values = c(
      SNP = "cadetblue3",
      MNP = "darkcyan",
      SV  = "slateblue4"
    )
  )
  ggsave(output, plot, width = 6, height = 4)
}