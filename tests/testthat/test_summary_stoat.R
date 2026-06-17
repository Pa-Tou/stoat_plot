library(testthat)

context("Summary STOAT function tests")

test_that("summary_stoat handles missing input file", {
  expect_error(summary_stoat("nonexistent.tsv"), "must be provided and exist")
})

test_that("summary_stoat validates parameters", {
  temp_file <- tempfile(fileext = ".tsv")
  file.create(temp_file)
  on.exit(unlink(temp_file))
  
  expect_error(summary_stoat(temp_file, number_top_var = -1), "must be positive numeric")
  expect_error(summary_stoat(temp_file, p_sig = 1.5), "must be >0 and <=1")
  expect_error(summary_stoat(temp_file, p_sig = 0), "must be >0 and <=1")
})
