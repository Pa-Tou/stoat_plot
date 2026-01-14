#' Dot Plot of Path Length Frequencies
#' @description Create a dot plot from a TSV file containing a `path_length` column.
#' 
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline
#' @importFrom utils read.table
#' 
#' @param input Path to the input TSV file.
#' @param output Path to save the output plot image.
#'
#' @return Saves a dot plot to the specified file.
#' @name path_length_distribution
#' @export

path_length_distribution <- function(input, output="paths_length_distribution.png") {

  # Read input
  df <- read.table(input, header = TRUE)

  if (!"TYPE" %in% colnames(df)) {
    stop("Column 'TYPE' not found in the data.")
  }

  # Split and flatten all TYPE values
  all_values <- unlist(sapply(as.character(df$TYPE), function(x) {
    vals <- unlist(strsplit(x, "[,/]+"))   # split on ',' & '/'
    as.numeric(vals)                       # convert all to numeric
  }))

  # Create data frame with counts
  freq_df <- as.data.frame(table(all_values))
  colnames(freq_df) <- c("Path_Length", "Frequency")
  freq_df$Path_Length <- as.numeric(as.character(freq_df$Path_Length))

  # ----------------- PLOT -----------------
  plot <- ggplot(freq_df, aes(x = Path_Length, y = Frequency)) +
    geom_point(size = 3, color = "steelblue") +
    labs(
      title = "Dot Plot of Path Length Frequency",
      x = "Path Length",
      y = "Frequency"
    ) +
    theme_bw() +
    theme(text = element_text(size = 14))

  # Save plot
  ggsave(output, plot, width = 6, height = 4)
}
