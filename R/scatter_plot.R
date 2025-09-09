#' Scatter Plot from Tab-Separated File
#'
#' @description
#' Generate a scatter plot from a tab-separated file, with support for optional
#' gzipped input, color grouping, and log transformation of the y-axis.
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
  out_file,
  title = "Title",
  x_label = "",
  y_label = "",
  x_col = 0,
  y_col = 1,
  color_col = -1,
  log_y = FALSE) {

  # Read data
  data <- utils::read.table(
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
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = group_label), alpha = 0.7) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label,
      color = if (!is.null(color_name)) paste0(color_name, ": count") else NULL
    ) +
    ggplot2::theme_minimal()

  # Apply log scale to y-axis if needed
  if (log_y) {
    p <- p + ggplot2::scale_y_log10()
  }

  # Save plot
  ggplot2::ggsave(out_file, plot = p, width = 12, height = 10, dpi = 400)
}
