
#' Dot Plot of snarl type histogram
#' @description Create a histogram plot from a TSV file containing a `path_length` column.
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
  plot <- ggplot2::ggplot(variant_counts, ggplot2::aes(x = Variant_Type, y = Count, fill = Variant_Type)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(title = "Snarl Type Distribution",
                  x = "Snarl Type",
                  y = "Count") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                   legend.position = "top")

  # Save plot
  ggplot2::ggsave(output, plot, width = 6, height = 4)
}