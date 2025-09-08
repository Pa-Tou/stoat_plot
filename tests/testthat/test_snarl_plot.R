test_that("snarl_type_histogram Binary creates PNG plot from valid input file", {
  # Locate test data file
  input_file <- system.file("data/snarl_paths", "binary_snarl_analyse.tsv", package = "StoatPlot")

  # Create temporary output file path
  output_file <- tempfile(fileext = ".png")

  # Run the function
  snarl_type_histogram(input = input_file, output = output_file)

  # Assert file is created
  expect_true(file.exists(output_file))

  # Optional: file is not empty
  expect_gt(file.info(output_file)$size, 0)
})

# test_that("snarl_type_histogram Quantitative creates PNG plot from valid input file", {
#   # Locate test data file
#   input_file <- system.file("data/snarl_paths", "quantitative_snarl_analyse.tsv", package = "StoatPlot")

#   # Create temporary output file path
#   output_file <- tempfile(fileext = ".png")

#   # Run the function
#   snarl_type_histogram(input = input_file, output = output_file)

#   # Assert file is created
#   expect_true(file.exists(output_file))

#   # Optional: file is not empty
#   expect_gt(file.info(output_file)$size, 0)
# })