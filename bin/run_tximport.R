#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tximport)
    library(optparse)
})

option_list <- list(
    make_option(
        "--tx2gene",
        type = "character"
    ),
    make_option(
        "--metadata",
        type = "character"
    ),
    make_option(
        "--counts",
        type = "character"
    ),
    make_option(
        "--txi",
        type = "character"
    ),
    make_option(
        "--summary",
        type = "character"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)

cat("====================================\n")
cat("RNAFlowX - TXIMPORT\n")
cat("====================================\n\n")


# ---------------------------------------------------------
# 1. Find Salmon quantification files
# ---------------------------------------------------------

quant_files <- list.files(
    path = ".",
    pattern = "^quant\\.sf$",
    recursive = TRUE,
    full.names = TRUE
)

if (length(quant_files) == 0) {
    stop(
        "No Salmon .sf files found in the working directory."
    )
}

cat(
    "Found",
    length(quant_files),
    "Salmon files:\n"
)

print(quant_files)


# ---------------------------------------------------------
# 2. Create sample names
# ---------------------------------------------------------

sample_names <- basename(dirname(quant_files))

names(quant_files) <- sample_names

cat("\nSample names:\n")
print(sample_names)


# ---------------------------------------------------------
# 3. Check expected samples
# ---------------------------------------------------------

metadata <- read.csv(
    opt$metadata,
    stringsAsFactors = FALSE
)

if (!"sample_id" %in% colnames(metadata)) {
    stop(
        "Metadata must contain a 'sample_id' column."
    )
}

expected_samples <- metadata$sample_id

missing_quant <- setdiff(
    expected_samples,
    sample_names
)

if (length(missing_quant) > 0) {

    stop(
        "Missing expected Salmon files: ",
        paste(
            missing_quant,
            collapse = ", "
        )
    )
}


# ---------------------------------------------------------
# 4. Read tx2gene
# ---------------------------------------------------------

cat("\nReading tx2gene...\n")

tx2gene <- read.delim(
    opt$tx2gene,
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t"
)

if (ncol(tx2gene) < 2) {
    stop(
        "tx2gene.tsv must contain at least two columns."
    )
}

tx2gene <- tx2gene[, 1:2]

colnames(tx2gene) <- c(
    "TXNAME",
    "GENEID"
)

cat(
    "tx2gene mappings:",
    nrow(tx2gene),
    "\n"
)


# ---------------------------------------------------------
# 5. Validate transcript compatibility
# ---------------------------------------------------------

cat("\nChecking transcript compatibility...\n")

example_quant <- quant_files[[1]]

quant_header <- read.delim(
    example_quant,
    header = TRUE,
    nrows = 5,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

quant_transcripts <- quant_header$Name

quant_transcripts_clean <- sub(
    "\\|.*$",
    "",
    quant_transcripts
)

matched <- sum(
    quant_transcripts_clean %in% tx2gene$TXNAME
)

cat(
    "Test transcripts:",
    length(quant_transcripts_clean),
    "\n"
)

cat(
    "Matched transcripts:",
    matched,
    "\n"
)

if (matched == 0) {

    stop(
        paste(
            "No transcript IDs from Salmon matched tx2gene.",
            "Check that the Salmon index and tx2gene file",
            "were generated from the same reference."
        )
    )
}


# ---------------------------------------------------------
# 6. Import Salmon quantifications
# ---------------------------------------------------------

cat("\nRunning tximport...\n")

txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreAfterBar = TRUE,
    countsFromAbundance = "no",
    dropInfReps = TRUE
)

cat(
    "\ntximport completed successfully.\n"
)

# ---------------------------------------------------------
# Save complete tximport object for DESeq2
# ---------------------------------------------------------

saveRDS(
    txi,
    file = opt$txi
)

cat(
    "Saved tximport object:",
    opt$txi,
    "\n"
)


# ---------------------------------------------------------
# 7. Create gene-level count matrix
# ---------------------------------------------------------

count_matrix <- as.data.frame(
    round(txi$counts)
)

if (nrow(count_matrix) == 0) {

    stop(
        "tximport produced an empty count matrix."
    )
}

write.csv(
    count_matrix,
    file = opt$counts,
    quote = FALSE
)

cat(
    "Count matrix dimensions:",
    nrow(count_matrix),
    "genes x",
    ncol(count_matrix),
    "samples\n"
)


# ---------------------------------------------------------
# 8. Read metadata
# ---------------------------------------------------------

cat("\nReading metadata...\n")

metadata <- read.csv(
    opt$metadata,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

write.csv(
    metadata,
    file = "metadata.csv",
    quote = FALSE
)

if (!"sample_id" %in% colnames(metadata)) {

    stop(
        "Metadata must contain a 'sample_id' column."
    )
}

if (!"condition" %in% colnames(metadata)) {

    stop(
        "Metadata must contain a 'condition' column."
    )
}

rownames(metadata) <- metadata$sample_id

metadata$sample_id <- NULL

cat("\nMetadata samples:\n")
print(rownames(metadata))


# ---------------------------------------------------------
# 9. Validate metadata against count matrix
# ---------------------------------------------------------

missing_metadata <- setdiff(
    colnames(count_matrix),
    rownames(metadata)
)

if (length(missing_metadata) > 0) {

    stop(
        "Count matrix samples missing from metadata: ",
        paste(
            missing_metadata,
            collapse = ", "
        )
    )
}

metadata <- metadata[
    colnames(count_matrix),
    ,
    drop = FALSE
]



# ---------------------------------------------------------
# 11. Create summary
# ---------------------------------------------------------

summary_lines <- c(
    "RNAFlowX TXIMPORT SUMMARY",
    "==========================",
    paste(
        "Samples:",
        ncol(count_matrix)
    ),
    paste(
        "Genes:",
        nrow(count_matrix)
    ),
    paste(
        "tx2gene mappings:",
        nrow(tx2gene)
    ),
    paste(
        "Counts method:",
        "Salmon NumReads"
    ),
    paste(
        "Inferential replicates:",
        "dropped"
    ),
    "",
    "Samples:",
    paste(
        colnames(count_matrix),
        collapse = ", "
    ),
    "",
    "Conditions:",
    paste(
        unique(metadata$condition),
        collapse = ", "
    )
)

writeLines(
    summary_lines,
    opt$summary
)


# ---------------------------------------------------------
# 12. Final message
# ---------------------------------------------------------

cat("\n====================================\n")
cat("TXIMPORT COMPLETED SUCCESSFULLY\n")
cat("====================================\n")

cat(
    "Count matrix:",
    opt$counts,
    "\n"
)

cat(
    "Summary:",
    opt$summary,
    "\n"
)
