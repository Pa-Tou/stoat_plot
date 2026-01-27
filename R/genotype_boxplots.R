#' Genotype boxplots for STOAT GWAS Results
#' @description Generates boxplots of phenotype by inferred genotype.
#'
#' @importFrom dplyr mutate filter select all_of left_join group_by summarise n na_if across
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_bw ggsave geom_abline
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @importFrom rlang .data
#' 
#' @param genotype_file Genotype file output of stoat vcf/graph
#' @param phenotype_file Path to the phenotype file use for the GWAS analysis.
#' @param node_start Node start boundary of the snarl [string]
#' @param node_end Snarl end boundary of the snarl [string]
#' @param output Path/Name to save the output plot image.
#'
#' @return Saves a genotype boxplots to the specified file.
#' @name genotype_boxplots
#' @export

genotype_boxplots <- function(genotype_file,
                              phenotype_file,
                              node_start,
                              node_end,
                              output = "boxplots.jpeg") {

  # -----------------------------
  # Read phenotype file
  # -----------------------------
  pheno_data <- read.table(
    phenotype_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "",
    check.names = FALSE
  )

  stopifnot(all(c("IID", "PHENO") %in% colnames(pheno_data)))

  # -----------------------------
  # Read STOAT genotype file (FORMAT-AWARE)
  # -----------------------------
  lines <- readLines(genotype_file)

  header_idx <- grep("^#START_NODE", lines)
  if (length(header_idx) == 0) {
    stop("START_NODE header not found in genotype file")
  }

  header <- sub("^#", "", lines[header_idx])
  col_names <- strsplit(header, "\t")[[1]]

  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]

  geno_data <- read.table(
    text = data_lines,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names,
    check.names = FALSE
  )

  # -----------------------------
  # Create snarl ID
  # -----------------------------
  geno_data <- geno_data %>%
    mutate(SNARL_ID = paste0(.data$START_NODE, .data$END_NODE))

  target_snarl <- paste0(node_start, node_end)

  geno_data <- geno_data %>%
    filter(.data$SNARL_ID == target_snarl)

  if (nrow(geno_data) == 0) {
    stop(sprintf("snarl_id %s not found in genotype file", target_snarl))
  }

  # -----------------------------
  # Identify genotype/sample columns
  # (all columns after SEQUENCES)
  # -----------------------------
  seq_idx <- which(colnames(geno_data) == "SEQUENCES")

  if (length(seq_idx) != 1) {
    stop("Column 'SEQUENCES' not found or not unique in genotype file.")
  }

  # all columns after SEQUENCES
  genotype_columns <- colnames(geno_data)[(seq_idx + 1):ncol(geno_data)]

  # Remove any derived columns that are not samples (e.g., SNARL_ID if exists)
  genotype_columns <- setdiff(genotype_columns, c("SNARL_ID"))

  if (length(genotype_columns) == 0) {
    stop("No genotype/sample columns detected after 'SEQUENCES'.")
  }

  # -----------------------------
  # Convert to long format
  # -----------------------------
  geno_long <- geno_data %>%
    mutate(across(all_of(genotype_columns), as.character)) %>%
    pivot_longer(
      cols = all_of(genotype_columns),
      names_to = "IID",
      values_to = "GT"
    ) %>%
    mutate(
      GT = na_if(.data$GT, "."),
      GT = as.numeric(.data$GT)
    ) %>%
    filter(!is.na(.data$GT))

  # N  checking R code for possible problems (7.1s)
  #   genotype_boxplots: no visible binding for global variable ‘GT’
  #   Undefined global functions or variables:
  #     GT

  # -----------------------------
  # Merge with phenotype
  # -----------------------------
  missing_in_pheno <- setdiff(genotype_columns, pheno_data$IID)
  if (length(missing_in_pheno) > 0) {
    stop(
      sprintf(
        "ERROR: Missing samples in phenotype file: %s",
        paste(missing_in_pheno, collapse = ", ")
      )
    )
  }

  missing_in_geno <- setdiff(pheno_data$IID, genotype_columns)
  if (length(missing_in_geno) > 0) {
    warning(
      sprintf(
        "WARNING: Phenotype samples not in genotype file: %s",
        paste(missing_in_geno, collapse = ", ")
      )
    )
  }

  merged_data <- geno_long %>%
    left_join(pheno_data, by = "IID") %>%
    filter(!is.na(.data$PHENO))

  # -----------------------------
  # Count genotypes
  # -----------------------------
  genotype_counts <- merged_data %>%
    group_by(.data$GT) %>%
    summarise(count = n(), .groups = "drop")

  merged_data <- merged_data %>%
    left_join(genotype_counts, by = "GT") %>%
    mutate(Genotype = paste0(.data$GT, "\n(", .data$count, ")"))

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(merged_data, aes(x = .data$Genotype, y = .data$PHENO)) +
    geom_violin(fill = "cadetblue3", alpha = 0.3) +
    geom_boxplot(
      width = 0.2,
      outlier.size = 2,
      outlier.colour = "red",
      alpha = 0.5,
      fill = "darkcyan"
    ) +
    labs(
      title = paste("Snarl:", target_snarl),
      x = "Genotype",
      y = "Phenotype"
    ) +
    theme_minimal(base_size = 14)

  # -----------------------------
  # Save
  # -----------------------------
  ggsave(
    filename = output,
    plot = p,
    device = "jpeg",
    width = 8,
    height = 6,
    dpi = 300
  )

  message("Saved plot: ", output)
}
