test_that("genotype_boxplots creates JPEG plots from phenotype and n_table files", {
  # Locate test phenotype file and genotype (n_table) directory
  phenotype_file <- system.file("data/phenotype", "path.tsv", package = "StoatPlot")
  n_table_dir <- system.file("data/regression", package = "StoatPlot")

  # Create temporary output directory
  output_dir <- tempfile("plots_")
  dir.create(output_dir)

  # Run the function
  genotype_boxplots(phenotype_file = phenotype_file, dir_path = n_table_dir, output = output_dir)

  # Check if any JPEG files were created
  plot_files <- list.files(output_dir, pattern = "\\.jpeg$", full.names = TRUE)

  # Expect at least one plot
  expect_true(length(plot_files) >= 1)

  # Optionally check file size is non-zero
  expect_true(all(file.info(plot_files)$size > 0))
})
