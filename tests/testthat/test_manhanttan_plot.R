test_that("manhattan_plot binary creates Manhattan plot image from input file", {

  input_file <- system.file("data/gwas", "binary_table_vcf.tsv", package = "StoatPlot")

  # Temporary output file
  output_file <- tempfile(fileext = "_binary.png")

  # Run the function
  manhattan_plot(input = input_file, output = output_file)

  # Check that PNG file was created and is not empty
  expect_true(file.exists(output_file))
  expect_gt(file.info(output_file)$size, 0)
})

# test_that("manhattan_plot quantitative creates Manhattan plot image from input file", {

#   input_file <- system.file("data/gwas", "quantitative_table_vcf.tsv", package = "StoatPlot")

#   # Temporary output file
#   output_file <- tempfile(fileext = "_quantitative.png")

#   # Run the function
#   manhattan_plot(input = input_file, output = output_file)

#   # Check that PNG file was created and is not empty
#   expect_true(file.exists(output_file))
#   expect_gt(file.info(output_file)$size, 0)
# })