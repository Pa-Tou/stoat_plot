#' Summary of Stoat snarl analysis
#'
#' @description Generate a summary of Stoat snarl output files.
#' @importFrom stats p.adjust median qchisq
#'
#' @param input Stoat snarl analysis file path (snarl gwas analyse) [string].
#' @param number_top_var Number of top variant print in the output top variant file (default: 100).
#' @param p_sig P-value threshold (default: 1e-5).
#' @param output Path/Name of the top variant output file (default: "top_variant.tsv").
#' @param apply_gc If TRUE, apply genomic control correction to p-values before reporting.
#'
#' @name summary_stoat
#' @export

summary_stoat <- function(input, number_top_var = 100, p_sig = 1e-5, output = "top_variant.tsv", apply_gc = FALSE) {
  if (is.null(input) || !file.exists(input)) stop("input must be provided and exist.")
  if (!is.numeric(number_top_var) || number_top_var <= 0) stop("number_top_var must be positive numeric.")
  if (!is.numeric(p_sig) || p_sig > 1 || p_sig <= 0) stop("p_sig must be >0 and <=1.")

  lines <- readLines(input)
  header_lines <- lines[grepl("^#", lines)]

  kv_headers <- header_lines[grepl(":", header_lines)]
  header_info <- if (length(kv_headers) > 0) setNames(sapply(strsplit(sub("^#", "", kv_headers), ":"), `[`, 2), sapply(strsplit(sub("^#", "", kv_headers), ":"), `[`, 1)) else list()

  ref_start <- which(header_lines == "#REFS")
  snarls_start <- which(header_lines == "#SNARLS")
  ref_names <- character(0)
  if (length(ref_start) == 1 && length(snarls_start) == 1 && snarls_start > ref_start) {
    ref_lines <- header_lines[(ref_start + 1):(snarls_start - 1)]
    ref_names <- trimws(sub("^#", "", ref_lines))
    ref_names <- ref_names[ref_names != ""]
  }

  data_lines <- lines[(which(lines == "#SNARLS") + 1):length(lines)]
  data_lines <- data_lines[!grepl("^#", data_lines)]
  col_names <- c("START_NODE", "END_NODE", "REF_INDEX", "START_OFFSET", "END_OFFSET", "DEPTH", "ALLELE_LENGTHS", "WALKS", "SEQUENCES", "P")
  df <- data.table::fread(paste(data_lines, collapse = "\n"), sep = "\t", header = FALSE, col.names = col_names, data.table = FALSE, quote = "", showProgress = FALSE)

  df$P <- suppressWarnings(as.numeric(df$P))
  total_variants <- nrow(df)
  p_significant <- sum(df$P < p_sig, na.rm = TRUE)

  n <- nrow(df)
  if (length(ref_names) > 0) {
    idx <- df$REF_INDEX + 1
    valid <- idx >= 1 & idx <= length(ref_names)
    df$CHR <- NA_character_
    df$CHR[valid] <- ref_names[idx[valid]]
    df$CHR[!valid] <- paste0("REF_", df$REF_INDEX[!valid])
  } else {
    df$CHR <- paste0("REF_", df$REF_INDEX)
  }
  variants_per_chr <- table(df$CHR)

  if ("P" %in% names(df) && length(df$P) == nrow(df)) df$P_BH <- stats::p.adjust(df$P, method = "BH") else df$P_BH <- rep(NA_real_, nrow(df))
  bh_significant_total <- sum(df$P_BH < p_sig, na.rm = TRUE)
  bh_significant_per_chr <- table(df$CHR[df$P_BH < p_sig])

  type_identification <- function(allele_lengths_str) {
    if (is.na(allele_lengths_str) || allele_lengths_str == ".") return(NA_character_)
    allele_lengths <- unlist(strsplit(allele_lengths_str, ",|/", perl = TRUE))
    allele_lengths <- as.integer(allele_lengths)
    if (length(allele_lengths) < 2) return(NA_character_)
    len_ref <- allele_lengths[1]
    alt_lengths <- allele_lengths[-1]
    variant_types <- character(0)
    for (len_alt in alt_lengths) {
      if (len_ref == 1 && len_alt == 1) variant_types <- c(variant_types, "SNP") else if (len_ref < 50 && len_alt < 50) variant_types <- c(variant_types, "MNP") else variant_types <- c(variant_types, "SV")
    }
    if ("SV" %in% variant_types) return("SV") else if ("MNP" %in% variant_types) return("MNP") else return("SNP")
  }

  df$VARIANT_TYPE <- sapply(df$ALLELE_LENGTHS, type_identification)
  variant_type_counts <- table(df$VARIANT_TYPE)
  variant_type_bh_significant <- table(df$VARIANT_TYPE[df$P_BH < p_sig])

  # genomic inflation factor
  chisq_stats <- stats::qchisq(1 - df$P, df = 1)
  lambda <- stats::median(chisq_stats, na.rm = TRUE) / stats::qchisq(0.5, df = 1)

  # Optionally apply genomic control and recompute BH
  if (apply_gc && !is.na(lambda) && lambda > 0) {
    chisq_gc <- chisq_stats / lambda
    df$P_gc <- pmax(1 - stats::pchisq(chisq_gc, df = 1), 1e-300)
    df$P_gc_BH <- stats::p.adjust(df$P_gc, method = "BH")
  }

  cat("====================================\n")
  cat("           Stoat Summary            \n")
  cat("====================================\n\n")

  if (length(kv_headers) > 0) {
    cat("Header Information:\n")
    for (nm in names(header_info)) cat(sprintf("  %s : %s\n", nm, header_info[[nm]]))
    cat("\n")
  }

  cat("Information Summary:\n")
  cat(sprintf("  Total variants                       : %d\n", total_variants))
  cat(sprintf("  Variants (with P < %.1e)           : %d\n", p_sig, p_significant))
  cat(sprintf("  Variants (with BH P < %.1e)        : %d\n", p_sig, bh_significant_total))
  cat(sprintf("  Genomic inflation factor (lambda)    : %.3f\n", lambda))
  cat("\n")

  if (!is.null(variants_per_chr)) {
    chr_names <- names(variants_per_chr)
    total_chr <- variants_per_chr
    p_sig_chr <- table(df$CHR[df$P < p_sig])
    p_sig_chr <- p_sig_chr[chr_names]
    p_sig_chr[is.na(p_sig_chr)] <- 0
    bh_sig_chr <- bh_significant_per_chr
    bh_sig_chr <- bh_sig_chr[chr_names]
    bh_sig_chr[is.na(bh_sig_chr)] <- 0
    chr_summary <- data.frame(Chromosome = chr_names, Total = as.integer(total_chr), P_lt_threshold = as.integer(p_sig_chr), BH_P_lt_threshold = as.integer(bh_sig_chr), check.names = FALSE, stringsAsFactors = FALSE)
    cat(sprintf("Per-Chromosome Variant Summary (Threshold: P < %.1e):\n", p_sig))
    print(chr_summary, row.names = FALSE)
    cat("\n")
  }

  if (!is.null(variant_type_counts)) {
    types <- names(variant_type_counts)
    total_type <- variant_type_counts[types]
    p_sig_type <- table(df$VARIANT_TYPE[df$P < p_sig])
    p_sig_type <- p_sig_type[types]
    p_sig_type[is.na(p_sig_type)] <- 0
    bh_sig_type <- variant_type_bh_significant[types]
    bh_sig_type[is.na(bh_sig_type)] <- 0
    type_summary <- data.frame(Variant_Type = types, Total = as.integer(total_type), P_lt_threshold = as.integer(p_sig_type), BH_P_lt_threshold = as.integer(bh_sig_type), check.names = FALSE, stringsAsFactors = FALSE)
    cat(sprintf("Variant Type Summary (Threshold: P < %.1e):\n", p_sig))
    print(type_summary, row.names = FALSE)
    cat("\n")
  }

  cat("====================================\n")

  df_top <- df[order(df$P, decreasing = FALSE), ]
  df_top <- head(df_top, number_top_var)
  cols_to_save <- c("START_NODE", "END_NODE", "REF_INDEX", "START_OFFSET", "END_OFFSET", "DEPTH", "ALLELE_LENGTHS", "SEQUENCES", "P")
  cols_to_save <- cols_to_save[cols_to_save %in% colnames(df_top)]
  write.table(df_top[, cols_to_save, drop = FALSE], file = output, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("The %d Top variants saved to: %s\n", number_top_var, output))
}
