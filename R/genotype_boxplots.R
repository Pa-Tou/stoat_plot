#' Genotype boxplots for STOAT GWAS Results
#' @description Generates boxplots of phenotype by inferred genotype.
#'
#' @param phenotype_file Path to the phenotype file use for the GWAS analysis.
#' @param dir_path Directory containing N-table files with genotype
#' @param output Path to save the output plot image.
#'
#' @return Saves a genotype boxplots to the specified file.
#' @name genotype_boxplots
#' @export
genotype_boxplots <- function(phenotype_file, dir_path, output="output_boxplots") {

  pheno_data <- read.table(phenotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  n_table_files <- list.files(path = dir_path, full.names = TRUE)

  for (n_table_file in n_table_files) {
    base_name <- tools::file_path_sans_ext(basename(n_table_file))

    # Read genotype table with custom header
    raw_header <- readLines(n_table_file, n = 1)
    header_cols <- strsplit(raw_header, "\t")[[1]]
    n_table <- read.table(n_table_file, header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
    colnames(n_table) <- header_cols
    colnames(n_table)[1] <- "IID"

    # Merge genotype with phenotype
    merged_data <- dplyr::left_join(n_table, pheno_data, by = "IID")

    genotype_columns <- setdiff(colnames(merged_data), c("IID", "PHENO"))

    # Inferred genotype from probabilities
    classify_genotype <- function(row) {
      probs <- suppressWarnings(as.numeric(row[genotype_columns]))
      if (all(is.na(probs))) return(NA)

      top_indices <- order(probs, decreasing = TRUE)[1:2]
      top_vals <- probs[top_indices]

      # Homozygous call: one allele has high probability
      if (!is.na(top_vals[1]) && top_vals[1] >= 0.95) {
        allele <- genotype_columns[top_indices[1]]
        return(paste0(allele, "/", allele))
      }

      # Heterozygous call: two alleles have balanced probabilities
      if (!is.na(top_vals[1]) && !is.na(top_vals[2]) &&
          top_vals[1] >= 0.4 && top_vals[2] >= 0.4) {
        allele1 <- genotype_columns[top_indices[1]]
        allele2 <- genotype_columns[top_indices[2]]
        return(paste0(allele1, "/", allele2))
      }

      return(NA)  # ambiguous case
    }

    merged_data$Genotype <- apply(merged_data, 1, classify_genotype)
    merged_data <- merged_data[!is.na(merged_data$PHENO) & !is.na(merged_data$Genotype), ]

    # Count genotypes
    genotype_counts <- dplyr::group_by(merged_data, Genotype) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop")

    # Append count to genotype label
    merged_data <- dplyr::left_join(merged_data, genotype_counts, by = "Genotype") %>%
      dplyr::mutate(Genotype = paste0(gsub("/", "\n", Genotype), "\n(", count, ")"))

    # Plot
    p <- ggplot2::ggplot(merged_data, ggplot2::aes(x = Genotype, y = PHENO)) +
      ggplot2::geom_violin(fill = "cadetblue3", alpha = 0.3) +
      ggplot2::geom_boxplot(width = 0.2, outlier.size = 2, outlier.colour = "red", alpha = 0.5, fill = "darkcyan") +
      ggplot2::labs(x = "Genotype", y = "Phenotype", title = paste("Boxplot -", base_name)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 14, color = "cadetblue4"),
        axis.title.y = ggplot2::element_text(size = 14, color = "cadetblue4"),
        plot.title = ggplot2::element_text(size = 16, color = "cadetblue4", face = "bold", hjust = 0.5)
      )

    # Save
    output_file <- file.path(output, paste0(base_name, "_boxplot.jpeg"))
    ggplot2::ggsave(output_file, plot = p, device = "jpeg", width = 8, height = 6, dpi = 300)
    message("Saved plot for ", base_name, " as ", output_file)
  }
}
