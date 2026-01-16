library(testthat)
library(StoatPlot)

test_that("StoatPlot main functions run without errors", {

  expect_error(summary_stoat("../../data/gwas/pg.gwas.tsv"), NA)
  expect_error(manhattan_plot("../../data/gwas/pg.gwas.tsv"), NA)
  expect_error(qq_plot("../../data/gwas/pg.gwas.tsv"), NA)
  expect_error(plot_pvalue_hist("../../data/gwas/pg.gwas.tsv", min = 1e-5, max = 0.5), NA)
  expect_error(genotype_boxplots(
    "../../data/genotype/pg.snarl.tsv",
    "../../data/phenotype/binary_phenotype_samples.tsv",
    "<4271", ">4260"), NA)
  expect_error(path_length_distribution("../../data/snarl_paths/binary_snarl_analyse.tsv"), NA)
  expect_error(snarl_type_histogram("../../data/snarl_paths/binary_snarl_analyse.tsv"), NA)
  expect_error(scatter_plot("../../data/snarl_paths/binary_snarl_analyse.tsv"), NA)

})