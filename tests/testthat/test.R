# devtools::test()

library(testthat)
library(StoatPlot)

test_that("StoatPlot main functions run without errors", {

  # Path to packaged test data
  test_data_dir <- system.file("extdata", package = "StoatPlot")
  expect_true(dir.exists(test_data_dir))

  # Temporary output directory (safe for CI)
  out_dir <- tempfile("stoatplot_test_")
  dir.create(out_dir)

  # Ensure cleanup even if test fails
  withr::defer(unlink(out_dir, recursive = TRUE), testthat::teardown_env())

  # -------------------------
  # Summary
  # -------------------------
  expect_error(
    summary_stoat(
      file.path(test_data_dir, "gwas", "pg.gwas.tsv"),
      output = file.path(out_dir, "summary.txt")
    ),
    NA
  )

  # -------------------------
  # Manhattan plot
  # -------------------------
  expect_error(
    manhattan_plot(
      file.path(test_data_dir, "gwas", "pg.gwas.tsv"),
      output = file.path(out_dir, "manhattan.png")
    ),
    NA
  )

  # -------------------------
  # Q-Q plot
  # -------------------------
  expect_error(
    qq_plot(
      file.path(test_data_dir, "gwas", "pg.gwas.tsv"),
      output = file.path(out_dir, "qq.png")
    ),
    NA
  )

  # -------------------------
  # P-value histogram
  # -------------------------
  expect_error(
    plot_pvalue_hist(
      file.path(test_data_dir, "gwas", "pg.gwas.tsv"),
      min = 1e-5,
      max = 0.5,
      output = file.path(out_dir, "pval_hist.png")
    ),
    NA
  )

  # -------------------------
  # Genotype boxplots
  # -------------------------
  expect_error(
    genotype_boxplots(
      file.path(test_data_dir, "genotype", "pg.snarl.tsv"),
      file.path(test_data_dir, "phenotype", "binary_phenotype_samples.tsv"),
      "<4271",
      ">4260",
      output = file.path(out_dir, "boxplots.png")
    ),
    NA
  )

  # -------------------------
  # Path length distribution
  # -------------------------
  expect_error(
    path_length_distribution(
      file.path(test_data_dir, "snarl_paths", "binary_snarl_analyse.tsv"),
      output = file.path(out_dir, "path_length.png")
    ),
    NA
  )

  # -------------------------
  # Snarl type histogram
  # -------------------------
  expect_error(
    snarl_type_histogram(
      file.path(test_data_dir, "snarl_paths", "binary_snarl_analyse.tsv"),
      output = file.path(out_dir, "snarl_type.png")
    ),
    NA
  )

  # -------------------------
  # Scatter plot
  # -------------------------
  expect_error(
    scatter_plot(
      file.path(test_data_dir, "snarl_paths", "binary_snarl_analyse.tsv"),
      output = file.path(out_dir, "scatter.png")
    ),
    NA
  )

})
