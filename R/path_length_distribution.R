#' Path Length Distribution Plots for STOAT
#' @description Generate multiple visualizations of path lengths from a TSV file with a `TYPE` column.
#'
#' @importFrom ggplot2 ggplot aes geom_point sec_axis labs theme_minimal theme element_text ggsave geom_line scale_y_continuous
#' @importFrom utils read.table
#' @importFrom stats ecdf
#'
#' @param input Path to the input TSV file.
#' @param min Minimum path length to include (default: 0).
#' @param max Maximum path length to include (default: Inf).
#' @param output Path to the output plot (default: "path_length_dot_ecdf.png").
#'
#' @return Saves four plots: dot plot, histogram, ECDF, and boxplot.
#' @name path_length_distribution_all
#' @export

path_length_distribution <- function(input, min = 0, max = Inf, output = "path_length_dot_ecdf.png") {

  # Read input
  df <- read.table(input, header = TRUE)
  if (!"TYPE" %in% colnames(df)) stop("Column 'TYPE' not found in data.")

  # Flatten path lengths
  all_values <- unlist(sapply(as.character(df$TYPE), function(x) {
    vals <- unlist(strsplit(x, "[,/]+"))
    as.numeric(vals)
  }))
  
  # Filter by min/max
  all_values <- all_values[!is.na(all_values) & all_values >= min & all_values <= max]
  if (length(all_values) == 0) stop("No path lengths within the specified range.")

  # Compute frequency table
  freq_df <- as.data.frame(table(all_values))
  colnames(freq_df) <- c("Path_Length", "Frequency")
  freq_df$Path_Length <- as.numeric(as.character(freq_df$Path_Length))
  
  # Compute ECDF
  ecdf_df <- data.frame(
    Path_Length = sort(unique(all_values)),
    Cumulative = ecdf(all_values)(sort(unique(all_values)))
  )

  # Scale ECDF to match frequency range for secondary axis
  max_freq <- max(freq_df$Frequency)
  ecdf_df$Cumulative_scaled <- ecdf_df$Cumulative * max_freq

  # ----------------- Combined Plot -----------------
  p <- ggplot() +
    # Dot plot for frequency
    geom_point(data = freq_df, aes(x = freq_df$Path_Length, y = freq_df$Frequency), 
              color = "cadetblue3", size = 3, alpha = 0.7) +
    # ECDF line scaled to frequency
    geom_line(data = ecdf_df, aes(x = freq_df$Path_Length, y = ecdf_df$Cumulative_scaled), 
              color = "darkcyan", linewidth = 1) +
    scale_y_continuous(
      name = "Frequency",
      sec.axis = sec_axis(~./max_freq, name = "Cumulative Proportion")
    ) +
    labs(
      title = paste0("Path Length Distribution & ECDF (", min, " - ", max, ")"),
      x = "Path Length"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title.y.left = element_text(color = "cadetblue3", face = "bold"),
      axis.title.y.right = element_text(color = "darkcyan", face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.text = element_text(color = "black", size = 12),
      panel.grid.minor = element_blank()
  )

  # Save plot
  ggsave(output, plot = p, width = 8, height = 5)

}
