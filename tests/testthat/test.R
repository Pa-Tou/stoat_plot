# devtools::test()

library(testthat)
library(StoatPlot)

test_that("StoatPlot main functions run without errors", {
  # Path to packaged test data
  test_data_dir <- system.file("extdata", package = "StoatPlot")
  expect_true(dir.exists(test_data_dir))

  # Temporary output directory (safe for CI)
  ## out_dir <- tempfile("stoatplot_test_")
  ## dir.create(out_dir)

  # Ensure cleanup even if test fails
  ## withr::defer(unlink(out_dir, recursive = TRUE), testthat::teardown_env())

  ## Summary
  expect_error(
    summary_stoat(file.path(test_data_dir, "stoat.binary.assoc.pvalues.tsv.gz")),
    NA)

  ## Manhattan plot
  expect_error(
    manhattan_plot(file.path(test_data_dir, "stoat.binary.assoc.pvalues.tsv.gz")),
    NA)

  ## Q-Q plot
  expect_error(
    qq_plot(file.path(test_data_dir, "stoat.binary.assoc.pvalues.tsv.gz")),
    NA)

  ## P-value histogram
  expect_error(
    plot_pvalue_hist(file.path(test_data_dir, "stoat.binary.assoc.pvalues.tsv.gz"),
                     min = 1e-5, max = 0.5, p_column="P_CHI2"),
    NA)

  ## Genotype boxplots
  expect_error(
    genotype_boxplots(genotype_file=file.path(test_data_dir, "snarl_genotypes.sorted.tsv.gz"),
                      phenotype=file.path(test_data_dir, "phenotype.quantitative.tsv"),
                      assoc=list(CHR="ref", START_OFFSET=53887,
                                 START_NODE='>1010', END_NODE='<1012')),
    NA)

  # -------------------------
  # Path length distribution
  # -------------------------
  ## expect_error(
  ##   path_length_distribution(
  ##     file.path(test_data_dir, "snarl_paths", "binary.snarl.tsv"),
  ##     output = file.path(out_dir, "path_length.png")
  ##   ),
  ##   NA
  ## )

  # -------------------------
  # Snarl type histogram
  # -------------------------
  ## expect_error(
  ##   snarl_type_histogram(
  ##     file.path(test_data_dir, "snarl_paths", "binary.snarl.tsv"),
  ##     output = file.path(out_dir, "snarl_type.png")
  ##   ),
  ##   NA
  ## )
})
