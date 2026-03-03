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
#' @param genotype_file Genotype path file output of stoat.
#' @param phenotype_file Phenotype path file use for the GWAS analysis.
#' @param node_start Start node boundary of the snarl.
#' @param node_end End node boundary of the snarl.
#' @param output Path/Name to save the output plot image.
#'
#' @return Saves a genotype boxplots to the specified file.
#' @name genotype_boxplots
#' @export

genotype_boxplots <- function(genotype_file,
                              phenotype_file,
                              node_start,
                              node_end,
                              output = "boxplots.png") {

  ## ---------------------------
  ## Input sanity checks
  ## ---------------------------
  if (!file.exists(genotype_file)) {
    stop("genotype_file does not exist: ", genotype_file)
  }

  if (!file.exists(phenotype_file)) {
    stop("phenotype_file does not exist: ", phenotype_file)
  }

  # -----------------------------
  # Read phenotype file
  # and validate format
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

  # stop is file is empty or only has header
  if (nrow(pheno_data) == 0) {
    stop("Phenotype file is empty or only has header.")
  }

  # -----------------------------
  # Read header Genotype file
  # and validate format
  # -----------------------------
  con <- file(genotype_file, open = "r")

  header_line <- NULL
  line_count <- 0

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break  # end of file
    line_count <- line_count + 1

    if (grepl("^#START_NODE\\b", line)) {
      header_line <- line
      break
    }
  }

  close(con)

  if (is.null(header_line)) {
    stop("header not found in genotype file.")
  }

  # Remove leading '#'
  header_line <- sub("^#", "", header_line)

  # Split columns
  col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]

  # Expected fixed header
  expected_cols <- c(
    "START_NODE",
    "END_NODE",
    "REF_INDEX",
    "START_OFFSET",
    "END_OFFSET",
    "DEPTH",
    "ALLELE_LENGTHS",
    "WALKS",
    "SEQUENCES"
  )

  # Check exact match
  if (!identical(col_names[1:9], expected_cols)) {
    stop(
      paste0(
        "Invalid genotype file format.\n",
        "Expected first 9 columns:\n",
        paste(expected_cols, collapse = "\t")
      )
    )
  }

  geno_data <- read.table(
    genotype_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = col_names,
    skip = line_count,
    comment.char = "",
    quote = "",
    fill = FALSE,
    check.names = FALSE
  )

  # -----------------------------
  # Create snarl ID
  # -----------------------------
  target_start <- node_start
  target_end   <- node_end

  # Filter directly (no SNARL_ID column needed)
  geno_data <- geno_data[
    geno_data$START_NODE == target_start &
    geno_data$END_NODE   == target_end,
    ,drop = FALSE
  ]

  if (nrow(geno_data) == 0) {
    stop(sprintf(
      "Snarl not found: %s -> %s. Try to swap start and end nodes.",
      target_start, target_end
    ))
  }

  if (nrow(geno_data) > 1) {
    stop(sprintf(
      "Duplicated Snarl found: %s -> %s. Please check the genotype file.",
      target_start, target_end
    ))
  }

  # -----------------------------
  # Identify genotype/sample columns
  # (all columns after SEQUENCES)
  # -----------------------------
  seq_idx <- match("SEQUENCES", names(geno_data))

  # all columns after SEQUENCES
  genotype_columns <- colnames(geno_data)[(seq_idx + 1):ncol(geno_data)]

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
      title = paste("Snarl:", target_start, target_end),
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
    width = 8,
    height = 6,
    dpi = 300
  )

  message("Saved plot: ", output)
}
