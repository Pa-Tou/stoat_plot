#' Genotype boxplots for STOAT GWAS Results
#' @description Generates boxplots of phenotype by inferred genotype.
#'
#' @param genotype_file Genotype file output of stoat vcf/graph
#' @param phenotype_file Path to the phenotype file use for the GWAS analysis.
#' @param snarl_id Snarl ID to plot
#' @param output Path/Name to save the output plot image.
#'
#' @return Saves a genotype boxplots to the specified file.
#' @name genotype_boxplots
#' @export

genotype_boxplots <- function(genotype_file,
                              phenotype_file,
                              snarl_id,
                              output = "boxplots.jpeg") {

  # -----------------------------
  # Read phenotype file
  # -----------------------------
  pheno_data <- read.table(
    phenotype_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  stopifnot(all(c("IID", "PHENO") %in% colnames(pheno_data)))

  # -----------------------------
  # Read STOAT genotype table
  # -----------------------------
  all_lines <- readLines(genotype_file)

  print(all_lines)

  # Find header line (WITH #)
  header_idx <- grep("^#START_NODE", all_lines)

  if (length(header_idx) == 0) {
    stop("START_NODE header not found in genotype file")
  }

  # Remove leading '#'
  header_line <- sub("^#", "", all_lines[header_idx])
  header_cols <- strsplit(header_line, "\t")[[1]]

  # Read data starting AFTER the header line
  geno_data <- read.table(
    genotype_file,
    header = FALSE,
    sep = "\t",
    skip = header_idx,
    stringsAsFactors = FALSE
  )

  colnames(geno_data) <- header_cols

  # -----------------------------
  # Create snarl ID (NO SEPARATOR)
  # -----------------------------
  geno_data <- geno_data %>%
    mutate(
      SNARL_ID = paste0(START_NODE, END_NODE)
    )

  # Filter selected snarl
  geno_data <- geno_data %>%
    filter(SNARL_ID == snarl_id)

  if (nrow(geno_data) == 0) {
    stop("snarl_id not found in genotype file")
  }

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
  geno_long <- geno_data %>%
    select(all_of(genotype_columns)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "IID",
      values_to = "GT"
    )

  # -----------------------------
  # Merge with phenotype
  # -----------------------------
  merged_data <- left_join(geno_long, pheno_data, by = "IID") %>%
    filter(!is.na(PHENO))

  # -----------------------------
  # Count genotypes
  # -----------------------------
  genotype_counts <- merged_data %>%
    group_by(GT) %>%
    summarise(count = n(), .groups = "drop")

  merged_data <- merged_data %>%
    left_join(genotype_counts, by = "GT") %>%
    mutate(Genotype = paste0(GT, "\n(", count, ")"))

  # -----------------------------
  # Plot
  # -----------------------------
  p <- ggplot(merged_data, aes(x = Genotype, y = PHENO)) +
    geom_violin(fill = "cadetblue3", alpha = 0.3) +
    geom_boxplot(
      width = 0.2,
      outlier.size = 2,
      outlier.colour = "red",
      alpha = 0.5,
      fill = "darkcyan"
    ) +
    labs(
      title = paste("Snarl:", snarl_id),
      x = "Genotype",
      y = "Phenotype"
    ) +
    theme_bw()

  # -----------------------------
  # Save output
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

