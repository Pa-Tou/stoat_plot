#' Summary of Stoat snarl analysis
#'
#' @description Generate a summary of Stoat snarl output files.
#' @importFrom utils read.table head write.table
#' @importFrom stats p.adjust median qchisq setNames
#'
#' @param input Stoat snarl analysis file path (snarl gwas analyse) [string].
#' @param number_top_var Number of top variant print in the output top variant file [string] (default : 100).
#' @param p_sig P-value threshold [numeric] (default : 1e-5).
#' @param output Path/Name of the top variant output file [string] (default : "top_variant.tsv").
#'
#' @name summary_stoat
#' @export

summary_stoat <- function(input, number_top_var=100, p_sig=1e-5, output="top_variant.tsv") {

    ## ---------------------------
    ## Input sanity checks
    ## ---------------------------
    if (is.null(input) || !file.exists(input)) stop("input must be provided and exist.")
    if (!is.numeric(number_top_var) || number_top_var <= 0) stop("number_top_var must be positive numeric.")
    if (!is.numeric(p_sig) || p_sig > 1 || p_sig <= 0) stop("p_sig must be >0 and <=1.")

    summary <- list()
    lines <- readLines(input)

    ## ---------------------------
    ## Extract header information
    ## ---------------------------
    header_lines <- lines[grepl("^#", lines)]

    # Key-value headers (e.g. #allele_size_limit:0)
    kv_headers <- header_lines[grepl(":", header_lines)]
    header_info <- setNames(
        sapply(strsplit(sub("^#", "", kv_headers), ":"), `[`, 2),
        sapply(strsplit(sub("^#", "", kv_headers), ":"), `[`, 1)
    )
    summary$header_info <- header_info

    ## ---------------------------
    ## Extract reference names
    ## ---------------------------
    ref_start <- which(header_lines == "#REFS")
    snarls_start <- which(header_lines == "#SNARLS")
    ref_names <- character(0)
    if (length(ref_start) == 1 && length(snarls_start) == 1 && snarls_start > ref_start) {
        ref_lines <- header_lines[(ref_start + 1):(snarls_start - 1)]
        ref_names <- trimws(sub("^#", "", ref_lines))
        ref_names <- ref_names[ref_names != ""]
    }

    ## ---------------------------
    ## Read variant table (force column names)
    ## ---------------------------
    data_lines <- lines[(which(lines == "#SNARLS")+1):length(lines)]
    data_lines <- data_lines[!grepl("^#", data_lines)]
    col_names <- c("START_NODE","END_NODE","REF_INDEX","START_OFFSET","END_OFFSET",
                    "DEPTH","ALLELE_LENGTHS","WALKS","SEQUENCES","P")
    df <- read.table(
        text = data_lines,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = col_names,
        quote = "",
        comment.char = ""
    )

    ## ---------------------------
    ## Ensure numeric P
    ## ---------------------------
    df$P <- suppressWarnings(as.numeric(df$P))

    ## ---------------------------
    ## Basic counts
    ## ---------------------------
    summary$total_variants <- nrow(df)
    summary$p_significant <- sum(df$P < p_sig, na.rm = TRUE)

    ## ---------------------------
    ## Per chromosome stats (SAFE)
    ## ---------------------------
    n <- nrow(df)
    df$CHR <- rep(NA_character_, n)
    if (length(ref_names) > 0) {
        idx <- df$REF_INDEX + 1
        valid <- idx >= 1 & idx <= length(ref_names)
        df$CHR[valid] <- ref_names[idx[valid]]
        df$CHR[!valid] <- paste0("REF_", df$REF_INDEX[!valid])
    } else {
        df$CHR <- paste0("REF_", df$REF_INDEX)
    }
    summary$variants_per_chr <- table(df$CHR)

    ## ---------------------------
    ## BH correction (SAFE)
    ## ---------------------------
    if ("P" %in% names(df) && length(df$P) == nrow(df)) {
        df$P_BH <- p.adjust(df$P, method = "BH")
    } else {
        df$P_BH <- rep(NA_real_, nrow(df))
    }
    summary$bh_significant_total <- sum(df$P_BH < p_sig, na.rm = TRUE)
    summary$bh_significant_per_chr <- table(df$CHR[df$P_BH < p_sig])

    # ---------------------------
    # Variant type identification based on ALLELE_LENGTHS
    # ---------------------------
    type_identification <- function(allele_lengths_str) {
        if (is.na(allele_lengths_str) || allele_lengths_str == ".") return(NA_character_)
        
        # Split multiple alleles; first part is ref, others are alt(s)
        # ALLELE_LENGTHS can be like: "1,1,5,5/6"
        allele_lengths <- unlist(strsplit(allele_lengths_str, ",|/", perl = TRUE))
        allele_lengths <- as.integer(allele_lengths)
        if (length(allele_lengths) < 2) return(NA_character_)
        
        len_ref <- allele_lengths[1]
        alt_lengths <- allele_lengths[-1]
        
        variant_types <- character(0)
        for (len_alt in alt_lengths) {
        if (len_ref == 1 && len_alt == 1) {
            variant_types <- c(variant_types, "SNP")
        } else if (len_ref < 50 && len_alt < 50) {
            variant_types <- c(variant_types, "MNP")
        } else {
            variant_types <- c(variant_types, "SV")
        }
        }
        
        # Priority order: SV > MNP > SNP
        if ("SV" %in% variant_types) {
        return("SV")
        } else if ("MNP" %in% variant_types) {
        return("MNP")
        } else {
        return("SNP")
        }
    }

    # ---------------------------
    # Apply to dataframe
    # ---------------------------
    df$VARIANT_TYPE <- sapply(df$ALLELE_LENGTHS, type_identification)

    summary$variant_type_counts <- table(df$VARIANT_TYPE)
    summary$variant_type_bh_significant <- table(df$VARIANT_TYPE[df$P_BH < p_sig])

    ## Store full processed table
    summary$variant_table <- df

    ## ---------------------------
    ## Genomic inflation factor λ
    ## ---------------------------
    chisq_stats <- qchisq(1 - df$P, df = 1)
    lambda <- median(chisq_stats, na.rm = TRUE) / qchisq(0.5, df = 1)
    summary$lambda <- lambda

    cat("====================================\n")
    cat("           Stoat Summary            \n")
    cat("====================================\n\n")

    ## ---------------------------
    ## Header information
    ## ---------------------------
    if (!is.null(summary$header_info) && length(summary$header_info) > 0) {
        cat("Header Information:\n")
        for (nm in names(summary$header_info)) {
        cat(sprintf("  %s : %s\n", nm, summary$header_info[[nm]]))
        }
        cat("\n")
    }

    ## ---------------------------
    ## Variant counts
    ## ---------------------------
    cat("Information Summary:\n")
    cat(sprintf("  Total variants                       : %d\n", summary$total_variants))
    cat(sprintf("  Variants (with P < %.1e)           : %d\n", p_sig, summary$p_significant))
    cat(sprintf("  Variants (with BH P < %.1e)        : %d\n", p_sig, summary$bh_significant_total))
    cat(sprintf("  Genomic inflation factor (lambda)    : %.3f\n", lambda))
    cat("\n")

    ## ---------------------------
    ## Per-chromosome summary table
    ## ---------------------------
    if (!is.null(summary$variants_per_chr) && !is.null(summary$variant_table)) {
        
        chr_names <- names(summary$variants_per_chr)
        
        # Total variants per chromosome
        total_chr <- summary$variants_per_chr
        
        # Variants with raw P < threshold
        p_sig_chr <- table(summary$variant_table$CHR[summary$variant_table$P < p_sig])
        # Ensure all chromosomes are present
        p_sig_chr <- p_sig_chr[chr_names]
        p_sig_chr[is.na(p_sig_chr)] <- 0
        
        # Variants with BH P < threshold
        bh_sig_chr <- summary$bh_significant_per_chr
        bh_sig_chr <- bh_sig_chr[chr_names]
        bh_sig_chr[is.na(bh_sig_chr)] <- 0
        
        # Combine into one table
        chr_summary <- data.frame(
        Chromosome = chr_names,
        Total = as.integer(total_chr),
        P_lt_threshold = as.integer(p_sig_chr),
        BH_P_lt_threshold = as.integer(bh_sig_chr),
        check.names = FALSE,
        stringsAsFactors = FALSE
        )
        
        cat(sprintf("Per-Chromosome Variant Summary (Threshold: P < %.1e):\n", p_sig))
        print(chr_summary, row.names = FALSE)
        cat("\n")
    }

    ## ---------------------------
    ## Variant type summary (all variants)
    ## ---------------------------
    if (!is.null(summary$variant_type_counts) && !is.null(summary$variant_table)) {
        
        types <- names(summary$variant_type_counts)
        
        # All variants
        total_type <- summary$variant_type_counts[types]
        
        # Variants with raw P < threshold
        p_sig_type <- table(summary$variant_table$VARIANT_TYPE[summary$variant_table$P < p_sig])
        p_sig_type <- p_sig_type[types]
        p_sig_type[is.na(p_sig_type)] <- 0
        
        # Variants with BH P < threshold
        bh_sig_type <- summary$variant_type_bh_significant[types]
        bh_sig_type[is.na(bh_sig_type)] <- 0
        
        # Combine into one table
        type_summary <- data.frame(
        Variant_Type = types,
        Total = as.integer(total_type),
        P_lt_threshold = as.integer(p_sig_type),
        BH_P_lt_threshold = as.integer(bh_sig_type),
        check.names = FALSE,
        stringsAsFactors = FALSE
        )
        
        cat(sprintf("Variant Type Summary (Threshold: P < %.1e):\n", p_sig))
        print(type_summary, row.names = FALSE)
        cat("\n")
    }

    ## ---------------------------
    ## Variant type per chromosome
    ## ---------------------------
    if (!is.null(summary$variant_table)) {
        chr_names <- unique(summary$variant_table$CHR)
        types <- c("SNP", "MNP", "SV")  # fixed order

        for (vt in types) {

        # Filter table for this variant type
        df_vt <- subset(summary$variant_table, summary$variant_table$VARIANT_TYPE == vt)
        if (nrow(df_vt) == 0) next  # skip if no variant of this type

        # Total per chromosome
        total_chr <- table(factor(df_vt$CHR, levels = chr_names))

        # Raw P < threshold
        p_chr <- table(factor(df_vt$CHR[df_vt$P < p_sig], levels = chr_names))

        # BH P < threshold
        bh_chr <- table(factor(df_vt$CHR[df_vt$P_BH < p_sig], levels = chr_names))

        # Combine into one data.frame
        vt_chr_summary <- data.frame(
            Chromosome = chr_names,
            Total = as.integer(total_chr),
            P_lt_threshold = as.integer(p_chr),
            BH_P_lt_threshold = as.integer(bh_chr),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        cat(sprintf("Variant Type: %s per Chromosome (Threshold: P < %.1e)\n", vt, p_sig))
        print(vt_chr_summary, row.names = FALSE)
        cat("\n")
        }
    }

    cat("====================================\n")

    # Order by P-value ascending
    df_top <- summary$variant_table[order(summary$variant_table$P, decreasing = FALSE), ]

    # Take top N variants
    df_top <- head(df_top, number_top_var)

    # Columns to save
    cols_to_save <- c("START_NODE", "END_NODE",
                    "REF_INDEX", "START_OFFSET", "END_OFFSET", 
                    "DEPTH", "ALLELE_LENGTHS", "SEQUENCES", "P")

    # Keep only columns that exist in df_top
    cols_to_save <- cols_to_save[cols_to_save %in% colnames(df_top)]

    # Save to file
    write.table(df_top[, cols_to_save, drop = FALSE], 
                file = output, sep = "\t",
                quote = FALSE, row.names = FALSE)
    cat(sprintf("The %d Top variants saved to: %s\n", number_top_var, output))
}
