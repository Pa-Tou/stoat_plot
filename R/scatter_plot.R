#' Scatter Plot from Tab-Separated File
#'
#' @description
#' Generate a scatter plot from a tab-separated file, with support for optional
#' gzipped input, color grouping, and log transformation of the y-axis.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_text ggsave geom_abline theme_minimal scale_y_log10
#' @importFrom utils read.table
#' 
#' @param input_file   Path to input file (.txt or .gz), tab-separated. Header is expected.
#' @param out_file     Path to save the output image (e.g., "output.png").
#' @param title        Plot title (default: "Title").
#' @param x_label      Label for the x-axis (default: column name from file).
#' @param y_label      Label for the y-axis (default: column name from file).
#' @param x_col        Index (0-based) of the x-axis column (default: 0).
#' @param y_col        Index (0-based) of the y-axis column (default: 1).
#' @param color_col    Index (0-based) of the column to group colors by (default: -1 = no color grouping).
#' @param log_y        Logical; whether to log-transform the y-axis (default: FALSE).
#'
#' @return
#' Saves a scatter plot image to the specified output path.
#'
#' @name scatter_plot
#' @export

scatter_plot <- function(
  input_file,
  out_file = "scatter_plot.png",
  title = "Title",
  x_label = "",
  y_label = "",
  x_col = 0,
  y_col = 1,
  color_col = -1,
  log_y = FALSE) {

  # Read data
  data <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )

  # Clean first column name if it starts with '#'
  colnames(data)[1] <- sub("^#", "", colnames(data)[1])

  # Get column names based on indices
  col_names <- colnames(data)
  x_name <- col_names[x_col + 1]
  y_name <- col_names[y_col + 1]
  color_name <- if (color_col != -1) col_names[color_col + 1] else NULL

  # Set axis labels if not provided
  if (x_label == "") x_label <- x_name
  if (y_label == "") y_label <- y_name

  # Extract relevant columns
  df <- data.frame(
    x = data[[x_name]],
    y = data[[y_name]],
    group = if (!is.null(color_name)) data[[color_name]] else "All",
    stringsAsFactors = FALSE
  )

  # Count per group for labeling
  group_counts <- table(df$group)
  df$group_label <- paste0(df$group, ": ", group_counts[df$group])

  # Build the plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(
      aes(color = group_label),
      alpha = 0.7,
      size = 2
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = if (!is.null(color_name)) paste0(color_name, ": count") else NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(
      values = rep(
        c("cadetblue3", "darkcyan"),
        length.out = length(unique(df$group_label))
      )
    )

  # Apply log scale to y-axis if needed
  if (log_y) {
    p <- p + scale_y_log10()
  }

  ggsave(out_file, plot = p, width = 12, height = 10, dpi = 400)
}
