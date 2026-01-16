
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

  df <- read.table(input, header = TRUE)

  # If TYPE is a comma-separated string per row, split and take max
  df$Variant_Type <- sapply(strsplit(as.character(df$TYPE), ","), function(path) {
    # Split by '/' and convert all to numeric
    values <- as.numeric(unlist(strsplit(path, "/")))
    values <- values[!is.na(values)]  # remove NA values

    if (length(values) == 0) {
      return(NA)  # no valid numbers in this path
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
  plot <- ggplot(variant_counts, aes(x = variant_counts$Variant_Type, y = variant_counts$Count, fill = variant_counts$Variant_Type)) +
    geom_bar(
      stat = "identity",
      position = "stack",
      alpha = 0.8
    ) +
    labs(
      title = "Snarl Type Distribution",
      x = "Snarl Type",
      y = "Count"
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
        SV  = "grey60"
      )
    )

  ggsave(output, plot, width = 6, height = 4)
}