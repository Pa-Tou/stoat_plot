library(testthat)
library(StoatPlot)

# Define the path to test data relative to the package root
test_data_dir <- "data"

test_that("StoatPlot main functions run without errors", {

  # Summary
  expect_error(
    summary_stoat(file.path(test_data_dir, "gwas", "pg.gwas.tsv")),
    NA
  )

  # Manhattan plot
  expect_error(
    manhattan_plot(file.path(test_data_dir, "gwas", "pg.gwas.tsv")),
    NA
  )

  # Q-Q plot
  expect_error(
    qq_plot(file.path(test_data_dir, "gwas", "pg.gwas.tsv")),
    NA
  )

  # P-value histogram
  expect_error(
    plot_pvalue_hist(file.path(test_data_dir, "gwas", "pg.gwas.tsv"),
                     min = 1e-5, max = 0.5),
    NA
  )

  # Genotype boxplots
  expect_error(
    genotype_boxplots(
      file.path(test_data_dir, "genotype", "pg.snarl.tsv"),
      file.path(test_data_dir, "phenotype", "binary_phenotype_samples.tsv"),
      "<4271", ">4260"),
    NA
  )

  # Path length distribution
  expect_error(
    path_length_distribution(file.path(test_data_dir, "snarl_paths", "binary_snarl_analyse.tsv")),
    NA
  )

  # Snarl type histogram
  expect_error(
    snarl_type_histogram(file.path(test_data_dir, "snarl_paths", "binary_snarl_analyse.tsv")),
    NA
  )

  # Scatter plot
  expect_error(
    scatter_plot(file.path(test_data_dir, "snarl_paths", "binary_snarl_analyse.tsv")),
    NA
  )

})
