#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(optparse)
})


# =========================================================
# Command-line arguments
# =========================================================

option_list <- list(

    make_option(
        "--counts",
        type = "character",
        help = "Input gene-level count matrix"
    ),

    make_option(
        "--metadata",
        type = "character",
        help = "Sample metadata CSV"
    ),

    make_option(
        "--outdir",
        type = "character",
        default = "deseq2_results",
        help = "Output directory"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)


# =========================================================
# Create output directory
# =========================================================

dir.create(
    opt$outdir,
    recursive = TRUE,
    showWarnings = FALSE
)


cat("====================================\n")
cat("RNAFlowX - DESeq2\n")
cat("====================================\n\n")


# =========================================================
# Read count matrix
# =========================================================

cat("Reading count matrix...\n")

counts <- read.csv(
    opt$counts,
    row.names = 1,
    check.names = FALSE
)

counts <- as.matrix(counts)

storage.mode(counts) <- "numeric"

cat(
    "Count matrix:",
    nrow(counts),
    "genes x",
    ncol(counts),
    "samples\n"
)


# =========================================================
# Read metadata
# =========================================================

cat("\nReading metadata...\n")

metadata <- read.csv(
    opt$metadata,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


# =========================================================
# Validate metadata
# =========================================================

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

print(
    rownames(metadata)
)


# =========================================================
# Validate sample matching
# =========================================================

count_samples <- colnames(counts)

metadata_samples <- rownames(metadata)


missing_metadata <- setdiff(
    count_samples,
    metadata_samples
)

extra_metadata <- setdiff(
    metadata_samples,
    count_samples
)


if (length(missing_metadata) > 0) {

    stop(
        "Samples present in count matrix but missing from metadata: ",
        paste(
            missing_metadata,
            collapse = ", "
        )
    )
}


if (length(extra_metadata) > 0) {

    warning(
        "Metadata contains extra samples: ",
        paste(
            extra_metadata,
            collapse = ", "
        )
    )
}


# Reorder metadata to match count matrix

metadata <- metadata[
    count_samples,
    ,
    drop = FALSE
]


cat("\nSample order verified:\n")

print(
    rownames(metadata)
)


# =========================================================
# Condition validation
# =========================================================

metadata$condition <- factor(
    metadata$condition
)

cat("\nConditions:\n")

print(
    table(metadata$condition)
)


if (nlevels(metadata$condition) != 2) {

    stop(
        "RNAFlowX currently expects exactly two conditions."
    )
}


# =========================================================
# Set control as reference
# =========================================================

if ("control" %in% levels(metadata$condition)) {

    metadata$condition <- relevel(
        metadata$condition,
        ref = "control"
    )

} else {

    warning(
        "Condition 'control' was not found. ",
        "Using the first condition as reference."
    )
}


conditions <- levels(
    metadata$condition
)

reference_condition <- conditions[1]
comparison_condition <- conditions[2]


cat(
    "\nReference condition:",
    reference_condition,
    "\n"
)

cat(
    "Comparison condition:",
    comparison_condition,
    "\n"
)


# =========================================================
# Create DESeq2 object
# =========================================================

cat("\nCreating DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = ~ condition
)


# =========================================================
# Filter low-count genes
# =========================================================

before_filter <- nrow(dds)

dds <- dds[
    rowSums(
        counts(dds) >= 10
    ) >= 2,
]

after_filter <- nrow(dds)

cat(
    "\nGenes before filtering:",
    before_filter,
    "\n"
)

cat(
    "Genes after filtering:",
    after_filter,
    "\n"
)


if (after_filter == 0) {

    stop(
        "No genes remain after count filtering."
    )
}


# =========================================================
# Run DESeq2
# =========================================================

cat("\nRunning DESeq2...\n")

dds <- DESeq(
    dds
)

cat(
    "DESeq2 completed successfully.\n"
)


# =========================================================
# Differential expression
# =========================================================

cat("\nCalculating differential expression...\n")

res <- results(
    dds,
    contrast = c(
        "condition",
        comparison_condition,
        reference_condition
    )
)

res <- as.data.frame(res)

res <- res[
    order(
        res$padj,
        na.last = NA
    ),
]


# =========================================================
# Save complete DE results
# =========================================================

write.csv(
    res,
    file = file.path(
        opt$outdir,
        "differential_expression.csv"
    ),
    quote = FALSE
)


# =========================================================
# Significant genes
# =========================================================

significant <- subset(
    res,
    !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)

write.csv(
    significant,
    file = file.path(
        opt$outdir,
        "significant_genes.csv"
    ),
    quote = FALSE
)


# =========================================================
# Normalized counts
# =========================================================

normalized_counts <- counts(
    dds,
    normalized = TRUE
)

write.csv(
    as.data.frame(normalized_counts),
    file = file.path(
        opt$outdir,
        "normalized_counts.csv"
    ),
    quote = FALSE
)


# =========================================================
# PCA
# =========================================================

cat("\nGenerating PCA plot...\n")

vsd <- vst(
    dds,
    blind = FALSE
)

pca_data <- plotPCA(
    vsd,
    intgroup = "condition",
    returnData = TRUE
)

percent_var <- round(
    100 * attr(
        pca_data,
        "percentVar"
    )
)

png(
    filename = file.path(
        opt$outdir,
        "PCA.png"
    ),
    width = 1800,
    height = 1600,
    res = 200
)

plot(
    pca_data$PC1,
    pca_data$PC2,
    pch = 19,
    xlab = paste0(
        "PC1: ",
        percent_var[1],
        "%"
    ),
    ylab = paste0(
        "PC2: ",
        percent_var[2],
        "%"
    ),
    main = "RNAFlowX PCA"
)

text(
    pca_data$PC1,
    pca_data$PC2,
    labels = rownames(pca_data),
    pos = 3,
    cex = 0.8
)

dev.off()


# =========================================================
# MA plot
# =========================================================

cat("Generating MA plot...\n")

png(
    filename = file.path(
        opt$outdir,
        "MA_plot.png"
    ),
    width = 1800,
    height = 1600,
    res = 200
)

plotMA(
    results(
        dds,
        contrast = c(
            "condition",
            comparison_condition,
            reference_condition
        )
    ),
    ylim = c(-5, 5),
    main = paste(
        comparison_condition,
        "vs",
        reference_condition
    )
)

dev.off()


# =========================================================
# Volcano plot
# =========================================================

cat("Generating volcano plot...\n")

volcano <- res

volcano$significance <- "Not significant"

volcano$significance[
    !is.na(volcano$padj) &
    volcano$padj < 0.05 &
    volcano$log2FoldChange >= 1
] <- "Upregulated"

volcano$significance[
    !is.na(volcano$padj) &
    volcano$padj < 0.05 &
    volcano$log2FoldChange <= -1
] <- "Downregulated"

volcano$neglog10padj <- -log10(
    volcano$padj
)

volcano$neglog10padj[
    is.infinite(volcano$neglog10padj)
] <- NA


png(
    filename = file.path(
        opt$outdir,
        "volcano_plot.png"
    ),
    width = 1800,
    height = 1600,
    res = 200
)

plot(
    volcano$log2FoldChange,
    volcano$neglog10padj,
    pch = 19,
    cex = 0.5,
    xlab = "log2 Fold Change",
    ylab = "-log10 adjusted p-value",
    main = paste(
        "Volcano:",
        comparison_condition,
        "vs",
        reference_condition
    )
)

abline(
    v = c(-1, 1),
    lty = 2
)

abline(
    h = -log10(0.05),
    lty = 2
)

dev.off()


# =========================================================
# Save DESeq2 object
# =========================================================

saveRDS(
    dds,
    file = file.path(
        opt$outdir,
        "dds.rds"
    )
)


# =========================================================
# Summary
# =========================================================

summary_lines <- c(

    "RNAFlowX DESeq2 SUMMARY",
    "========================",

    paste(
        "Genes before filtering:",
        before_filter
    ),

    paste(
        "Genes after filtering:",
        after_filter
    ),

    paste(
        "Samples:",
        ncol(counts)
    ),

    paste(
        "Reference:",
        reference_condition
    ),

    paste(
        "Comparison:",
        comparison_condition
    ),

    paste(
        "Significant genes:",
        nrow(significant)
    ),

    "",

    "Samples:",

    paste(
        colnames(counts),
        collapse = ", "
    ),

    "",

    "Conditions:",

    paste(
        unique(
            as.character(
                metadata$condition
            )
        ),
        collapse = ", "
    )
)

writeLines(
    summary_lines,
    file.path(
        opt$outdir,
        "deseq2_summary.txt"
    )
)


# =========================================================
# Final message
# =========================================================

cat("\n====================================\n")
cat("DESEQ2 COMPLETED SUCCESSFULLY\n")
cat("====================================\n")

cat(
    "Reference:",
    reference_condition,
    "\n"
)

cat(
    "Comparison:",
    comparison_condition,
    "\n"
)

cat(
    "Significant genes:",
    nrow(significant),
    "\n"
)

cat(
    "Output directory:",
    opt$outdir,
    "\n"
)