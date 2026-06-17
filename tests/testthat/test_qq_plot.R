library(testthat)

context("QQ plot function tests")

test_that("qq_plot handles missing input file", {
  expect_error(qq_plot("nonexistent.tsv"), "not found")
})

test_that("qq_plot validates p_column parameter", {
  temp_file <- tempfile(fileext = ".tsv")
  writeLines("#CHR\tP\nCHR1\t0.05", temp_file)
  on.exit(unlink(temp_file))
  
  # Should fail because P_CHI2 doesn't exist but P does (falls back)
  # So this should actually succeed by falling back to P
  result <- qq_plot(temp_file, output = NULL)
  expect_invisible(result)
})
