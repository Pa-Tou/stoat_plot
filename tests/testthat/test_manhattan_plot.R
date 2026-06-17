library(testthat)

context("Manhattan plot function tests")

test_that("manhattan_plot handles missing input file", {
  expect_error(manhattan_plot("nonexistent.tsv"), "not found")
})

test_that("manhattan_plot returns invisibly when output is NULL", {
  # Create minimal mock TSV with required header
  temp_file <- tempfile(fileext = ".tsv")
  writeLines("#CHR\tSTART_OFFSET\tP", temp_file)
  on.exit(unlink(temp_file))
  
  result <- manhattan_plot(temp_file, output = NULL)
  expect_invisible(result)
})
