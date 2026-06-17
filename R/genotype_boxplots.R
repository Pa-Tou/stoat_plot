#' Genotype boxplots for STOAT GWAS Results
#' @description Generates boxplots of phenotype by inferred genotype.
#'
#' @importFrom dplyr mutate filter left_join group_by summarise n across
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_minimal ggsave
#' @importFrom rlang .data
#'
#' @param genotype_file Genotype file output of stoat vcf/graph
#' @param phenotype_file Path to the phenotype file used for the GWAS analysis.
#' @param node_start Node start boundary of the snarl [string]
#' @param node_end Snarl end boundary of the snarl [string]
#' @param output Path/Name to save the output plot image.
#' @param show_test If TRUE, compute and display a simple ANOVA p-value comparing genotypes.
#'
#' @return Saves a genotype boxplot to the specified file and returns the ggplot object invisibly.
#' @name genotype_boxplots
#' @export

genotype_boxplots <- function(genotype_file,
                              phenotype_file,
                              node_start,
                              node_end,
                              output = "boxplots.jpeg",
                              show_test = TRUE) {
  pheno_data <- data.table::fread(phenotype_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  stopifnot(all(c("IID", "PHENO") %in% colnames(pheno_data)))

  lines <- readLines(genotype_file)
  header_idx <- grep("^#START_NODE", lines)
  if (length(header_idx) == 0) stop("START_NODE header not found in genotype file")
  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  geno <- data.table::fread(paste(data_lines, collapse = "\n"), sep = "\t", header = FALSE, col.names = col_names, data.table = FALSE, check.names = FALSE)

  geno <- dplyr::mutate(geno, SNARL_ID = paste0(.data$START_NODE, .data$END_NODE))
  target_snarl <- paste0(node_start, node_end)
  geno <- dplyr::filter(geno, .data$SNARL_ID == target_snarl)
  if (nrow(geno) == 0) stop(sprintf("snarl_id %s not found in genotype file", target_snarl))

  seq_idx <- which(colnames(geno) == "SEQUENCES")
  if (length(seq_idx) != 1) stop("Column 'SEQUENCES' not found or not unique in genotype file.")
  genotype_columns <- colnames(geno)[(seq_idx + 1):ncol(geno)]
  genotype_columns <- setdiff(genotype_columns, c("SNARL_ID"))
  if (length(genotype_columns) == 0) stop("No genotype/sample columns detected after 'SEQUENCES'.")

  geno_long <- geno %>%
    dplyr::mutate(dplyr::across(all_of(genotype_columns), as.character)) %>%
    tidyr::pivot_longer(cols = all_of(genotype_columns), names_to = "IID", values_to = "GT") %>%
    dplyr::mutate(GT = dplyr::na_if(.data$GT, "."), GT = as.numeric(.data$GT)) %>%
    dplyr::filter(!is.na(.data$GT))

  missing_in_pheno <- setdiff(genotype_columns, pheno_data$IID)
  if (length(missing_in_pheno) > 0) stop(sprintf("ERROR: Missing samples in phenotype file: %s", paste(missing_in_pheno, collapse = ", ")))
  missing_in_geno <- setdiff(pheno_data$IID, genotype_columns)
  if (length(missing_in_geno) > 0) warning(sprintf("WARNING: Phenotype samples not in genotype file: %s", paste(missing_in_geno, collapse = ", ")))

  merged_data <- dplyr::left_join(geno_long, pheno_data, by = "IID") %>% dplyr::filter(!is.na(.data$PHENO))

  genotype_counts <- merged_data %>% dplyr::group_by(.data$GT) %>% dplyr::summarise(count = n(), .groups = "drop")
  merged_data <- merged_data %>% dplyr::left_join(genotype_counts, by = "GT") %>% dplyr::mutate(Genotype = paste0(.data$GT, "\n(", .data$count, ")"))

  p <- ggplot2::ggplot(merged_data, ggplot2::aes(x = .data$Genotype, y = .data$PHENO)) +
    ggplot2::geom_violin(fill = "cadetblue3", alpha = 0.3) +
    ggplot2::geom_boxplot(width = 0.2, outlier.size = 2, outlier.colour = "red", alpha = 0.5, fill = "darkcyan") +
    ggplot2::labs(title = paste("Snarl:", target_snarl), x = "Genotype", y = "Phenotype") +
    ggplot2::theme_minimal(base_size = 14)

  # optional simple test
  if (show_test) {
    try({
      anov <- stats::aov(PHENO ~ factor(GT), data = merged_data)
      pval <- summary(anov)[[1]]["Pr(>F)"][1]
      p <- p + ggplot2::labs(subtitle = paste0("ANOVA p = ", signif(pval, 3)))
    }, silent = TRUE)
  }

  ggplot2::ggsave(filename = output, plot = p, device = "jpeg", width = 8, height = 6, dpi = 300)
  message("Saved plot: ", output)
  invisible(p)
}
