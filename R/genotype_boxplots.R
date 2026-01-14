#' Genotype boxplots for STOAT GWAS Results
#' @description Generates boxplots of phenotype by inferred genotype.
#'
#' @param genotype_file Genotype file output of stoat vcf/graph
#' @param phenotype_file Path to the phenotype file use for the GWAS analysis.
#' @param node_start Node start boundary of the snarl
#' @param node_end Snarl end boundary of the snarl

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
    comment.char = ""
  )

  stopifnot(all(c("IID", "PHENO") %in% colnames(pheno_data)))

  # -----------------------------
  # Read STOAT genotype table
  # -----------------------------
  all_lines <- readLines(genotype_file)

  header_idx <- grep("^#START_NODE", all_lines)
  if (length(header_idx) == 0) {
    stop("START_NODE header not found in genotype file")
  }

  header_line <- sub("^#", "", all_lines[header_idx])
  header_cols <- strsplit(header_line, "\t")[[1]]

  geno_data <- read.table(
    genotype_file,
    header = FALSE,
    sep = "\t",
    skip = header_idx,
    stringsAsFactors = FALSE
  )

  colnames(geno_data) <- header_cols

  # -----------------------------
  # Create snarl ID
  # -----------------------------
  geno_data <- geno_data |>
    dplyr::mutate(
      SNARL_ID = paste0(
        gsub("[<>]", "", START_NODE),
        "_",
        gsub("[<>]", "", END_NODE)
      )
    )

  target_snarl <- paste0(node_start, "_", node_end)

  filtered <- geno_data |>
    dplyr::filter(SNARL_ID == target_snarl)

  if (nrow(filtered) == 0) {
    message(sprintf(
      "snarl_id %s not found, trying reversed order",
      target_snarl
    ))

    target_snarl_rev <- paste0(node_end, "_", node_start)

    filtered <- geno_data |>
      dplyr::filter(SNARL_ID == target_snarl_rev)

    if (nrow(filtered) == 0) {
      stop(sprintf(
        "snarl_id %s or %s not found in genotype file",
        target_snarl,
        target_snarl_rev
      ))
    }
  }

  geno_data <- filtered

  # -----------------------------
  # Identify genotype columns
  # -----------------------------
  fixed_cols <- c(
    "START_NODE", "END_NODE", "REF", "START_OFFSET", "END_OFFSET",
    "DEPTH", "ALLELE_LENGTHS", "WALKS", "SEQUENCES", "SNARL_ID"
  )

  genotype_columns <- setdiff(colnames(geno_data), fixed_cols)

  # -----------------------------
  # Convert to sample-long format
  # -----------------------------
  geno_long <- geno_data |>
    dplyr::select(dplyr::all_of(genotype_columns)) |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "IID",
      values_to = "GT"
    )

  # -----------------------------
  # Merge genotype with phenotype
  # -----------------------------
  missing_in_pheno <- setdiff(genotype_columns, pheno_data$IID)
  if (length(missing_in_pheno) > 0) {
    stop(sprintf(
      "ERROR: The following genotype samples are missing in the phenotype file: %s",
      paste(missing_in_pheno, collapse = ", ")
    ))
  }

  missing_in_geno <- setdiff(pheno_data$IID, genotype_columns)
  if (length(missing_in_geno) > 0) {
    warning(sprintf(
      "WARNING: The following phenotype samples are not present in the genotype file: %s",
      paste(missing_in_geno, collapse = ", ")
    ))
  }

  merged_data <- geno_long |>
    dplyr::left_join(pheno_data, by = "IID") |>
    dplyr::filter(!is.na(PHENO))

  # -----------------------------
  # Count genotypes
  # -----------------------------
  genotype_counts <- merged_data |>
    dplyr::group_by(GT) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  merged_data <- merged_data |>
    dplyr::left_join(genotype_counts, by = "GT") |>
    dplyr::mutate(Genotype = paste0(GT, "\n(", count, ")"))

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot2::ggplot(
    merged_data,
    ggplot2::aes(x = Genotype, y = PHENO)
  ) +
    ggplot2::geom_violin(fill = "cadetblue3", alpha = 0.3) +
    ggplot2::geom_boxplot(
      width = 0.2,
      outlier.size = 2,
      outlier.colour = "red",
      alpha = 0.5,
      fill = "darkcyan"
    ) +
    ggplot2::labs(
      title = paste("Snarl:", target_snarl),
      x = "Genotype",
      y = "Phenotype"
    ) +
    ggplot2::theme_bw()

  # -----------------------------
  # Save output
  # -----------------------------
  ggplot2::ggsave(
    filename = output,
    plot = p,
    device = "jpeg",
    width = 8,
    height = 6,
    dpi = 300
  )

  message("Saved plot: ", output)
}
