library(testthat)

context("Basic StoatPlot smoke tests")

test_that("functions error on missing files", {
  expect_error(manhattan_plot("nonexistent.tsv"), "not found")
  expect_error(qq_plot("nonexistent.tsv"), "not found")
  expect_error(genotype_boxplots("nonexistent.tsv", "nonexistent2.tsv", "1", "2"), "not found")
  expect_error(summary_stoat("nonexistent.tsv"), "must be provided and exist")
})
