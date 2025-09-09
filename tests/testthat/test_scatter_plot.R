test_that("scatter_plot creates a PNG scatter plot from a tab-separated file", {
  # Locate test input file (supports .txt or .gz)
  input_file <- system.file("extdata", "scatter_data.tsv", package = "StoatPlot")
  expect_true(file.exists(input_file))

  # Create temporary output file
  output_file <- tempfile(fileext = ".png")

  # Run the function
  StoatPlot::scatter_plot(
    input_file = input_file,
    out_file = output_file,
    title = "Test Scatter",
    x_label = "X Axis",
    y_label = "Y Axis",
    x_col = 0,
    y_col = 1,
    color_col = -1,
    log_y = FALSE
  )

  # Check if PNG file was created
  expect_true(file.exists(output_file))

  # Check that file is non-empty
  expect_true(file.info(output_file)$size > 0)
})
