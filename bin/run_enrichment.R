#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop(
        "Usage: Rscript run_enrichment.R ",
        "<differential_expression.csv> <outdir>"
    )
}

input_file <- args[1]
outdir <- args[2]

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("====================================\n")
cat("RNAFlowX - Functional Enrichment\n")
cat("====================================\n\n")

# =========================================================
# Read DESeq2 results
# =========================================================

cat("Reading DESeq2 results...\n")

res <- read.csv(
    input_file,
    row.names = 1,
    check.names = FALSE
)

res <- res[
    !is.na(res$padj) &
    !is.na(res$log2FoldChange) &
    !is.na(res$stat),
]

cat(
    "Genes available:",
    nrow(res),
    "\n"
)

# =========================================================
# Remove Ensembl version suffix
# =========================================================

res$ensembl_id <- sub(
    "\\..*$",
    "",
    rownames(res)
)

res <- res[
    !duplicated(res$ensembl_id),
]

cat(
    "Unique Ensembl genes:",
    nrow(res),
    "\n"
)

# =========================================================
# Convert Ensembl -> Entrez
# =========================================================

cat("\nMapping Ensembl IDs to Entrez IDs...\n")

gene_map <- bitr(
    res$ensembl_id,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

gene_map <- gene_map[
    !duplicated(gene_map$ENSEMBL),
]

cat(
    "Mapped genes:",
    nrow(gene_map),
    "\n"
)

# =========================================================
# Add Entrez IDs
# =========================================================

res$ENTREZID <- gene_map$ENTREZID[
    match(
        res$ensembl_id,
        gene_map$ENSEMBL
    )
]

mapped_res <- res[
    !is.na(res$ENTREZID),
]

cat(
    "Genes retained after mapping:",
    nrow(mapped_res),
    "\n"
)

# =========================================================
# Define significant genes
# =========================================================

significant <- mapped_res[
    mapped_res$padj < 0.05 &
    abs(mapped_res$log2FoldChange) >= 1,
]

upregulated <- significant[
    significant$log2FoldChange >= 1,
]

downregulated <- significant[
    significant$log2FoldChange <= -1,
]

cat(
    "\nSignificant genes:",
    nrow(significant),
    "\n"
)

cat(
    "Upregulated:",
    nrow(upregulated),
    "\n"
)

cat(
    "Downregulated:",
    nrow(downregulated),
    "\n"
)

# =========================================================
# Output directories
# =========================================================

go_bp_dir <- file.path(
    outdir,
    "GO",
    "BP"
)

go_mf_dir <- file.path(
    outdir,
    "GO",
    "MF"
)

go_cc_dir <- file.path(
    outdir,
    "GO",
    "CC"
)

kegg_dir <- file.path(
    outdir,
    "KEGG"
)

plots_dir <- file.path(
    outdir,
    "plots"
)

dirs <- c(
    go_bp_dir,
    go_mf_dir,
    go_cc_dir,
    kegg_dir,
    plots_dir
)

for (d in dirs) {
    dir.create(
        d,
        recursive = TRUE,
        showWarnings = FALSE
    )
}

# =========================================================
# Helper function for ORA
# =========================================================

run_go <- function(
    genes,
    ontology,
    output_file
) {

    if (length(genes) == 0) {
        return(NULL)
    }

    result <- enrichGO(
        gene = genes,
        universe = mapped_res$ENTREZID,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )

    if (!is.null(result)) {

        write.csv(
            as.data.frame(result),
            output_file,
            row.names = FALSE
        )
    }

    result
}

# =========================================================
# GO - All significant genes
# =========================================================

cat("\nRunning GO enrichment - all significant genes...\n")

go_bp_all <- run_go(
    significant$ENTREZID,
    "BP",
    file.path(
        go_bp_dir,
        "all_significant.csv"
    )
)

go_mf_all <- run_go(
    significant$ENTREZID,
    "MF",
    file.path(
        go_mf_dir,
        "all_significant.csv"
    )
)

go_cc_all <- run_go(
    significant$ENTREZID,
    "CC",
    file.path(
        go_cc_dir,
        "all_significant.csv"
    )
)

# =========================================================
# GO - Upregulated
# =========================================================

cat("Running GO enrichment - upregulated genes...\n")

go_bp_up <- run_go(
    upregulated$ENTREZID,
    "BP",
    file.path(
        go_bp_dir,
        "upregulated.csv"
    )
)

go_mf_up <- run_go(
    upregulated$ENTREZID,
    "MF",
    file.path(
        go_mf_dir,
        "upregulated.csv"
    )
)

go_cc_up <- run_go(
    upregulated$ENTREZID,
    "CC",
    file.path(
        go_cc_dir,
        "upregulated.csv"
    )
)

# =========================================================
# GO - Downregulated
# =========================================================

cat("Running GO enrichment - downregulated genes...\n")

go_bp_down <- run_go(
    downregulated$ENTREZID,
    "BP",
    file.path(
        go_bp_dir,
        "downregulated.csv"
    )
)

go_mf_down <- run_go(
    downregulated$ENTREZID,
    "MF",
    file.path(
        go_mf_dir,
        "downregulated.csv"
    )
)

go_cc_down <- run_go(
    downregulated$ENTREZID,
    "CC",
    file.path(
        go_cc_dir,
        "downregulated.csv"
    )
)

# =========================================================
# KEGG
# =========================================================

cat("\nRunning KEGG enrichment...\n")

kegg_all <- enrichKEGG(
    gene = significant$ENTREZID,
    universe = mapped_res$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

write.csv(
    as.data.frame(kegg_all),
    file.path(
        kegg_dir,
        "all_significant.csv"
    ),
    row.names = FALSE
)

kegg_up <- enrichKEGG(
    gene = upregulated$ENTREZID,
    universe = mapped_res$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

write.csv(
    as.data.frame(kegg_up),
    file.path(
        kegg_dir,
        "upregulated.csv"
    ),
    row.names = FALSE
)

kegg_down <- enrichKEGG(
    gene = downregulated$ENTREZID,
    universe = mapped_res$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

write.csv(
    as.data.frame(kegg_down),
    file.path(
        kegg_dir,
        "downregulated.csv"
    ),
    row.names = FALSE
)

# =========================================================
# Save mapped results
# =========================================================

dir.create(
    file.path(
        outdir,
        "annotation"
    ),
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    mapped_res,
    file.path(
        outdir,
        "annotation",
        "mapped_DESeq2_results.csv"
    ),
    row.names = TRUE
)

write.csv(
    gene_map,
    file.path(
        outdir,
        "annotation",
        "ensembl_to_entrez.csv"
    ),
    row.names = FALSE
)

# =========================================================
# GSEA
# =========================================================

cat("\nPreparing GSEA ranking...\n")

gene_list <- mapped_res$stat

names(gene_list) <- mapped_res$ENTREZID

gene_list <- gene_list[
    !is.na(gene_list)
]

gene_list <- sort(
    gene_list,
    decreasing = TRUE
)

gene_list <- gene_list[
    !duplicated(names(gene_list))
]

cat(
    "Genes in GSEA ranking:",
    length(gene_list),
    "\n"
)

gsea_go <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
)

dir.create(
    file.path(
        outdir,
        "GSEA"
    ),
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    as.data.frame(gsea_go),
    file.path(
        outdir,
        "GSEA",
        "GO_BP_GSEA.csv"
    ),
    row.names = FALSE
)

# =========================================================
# Visualization
# =========================================================

cat("\nGenerating enrichment plots...\n")

plot_enrichment <- function(
    result,
    output_file,
    title
) {

    if (is.null(result)) {
        cat("Skipping:", title, "- no enrichment result\n")
        return(invisible(NULL))
    }

    result_df <- as.data.frame(result)

    if (nrow(result_df) == 0) {
        cat("Skipping:", title, "- no enriched terms\n")
        return(invisible(NULL))
    }

    p <- dotplot(
        result,
        showCategory = 15
    ) +
        ggtitle(title)

    ggsave(
        filename = output_file,
        plot = p,
        width = 10,
        height = 7,
        dpi = 300
    )

    cat("Created:", output_file, "\n")
}

plot_enrichment(
    go_bp_all,
    file.path(
        plots_dir,
        "GO_BP_all_significant_dotplot.png"
    ),
    "GO Biological Process - All Significant Genes"
)

plot_enrichment(
    go_bp_up,
    file.path(
        plots_dir,
        "GO_BP_upregulated_dotplot.png"
    ),
    "GO Biological Process - Upregulated Genes"
)

plot_enrichment(
    go_bp_down,
    file.path(
        plots_dir,
        "GO_BP_downregulated_dotplot.png"
    ),
    "GO Biological Process - Downregulated Genes"
)

plot_enrichment(
    kegg_all,
    file.path(
        plots_dir,
        "KEGG_all_significant_dotplot.png"
    ),
    "KEGG Pathway Enrichment - All Significant Genes"
)

plot_enrichment(
    kegg_up,
    file.path(
        plots_dir,
        "KEGG_upregulated_dotplot.png"
    ),
    "KEGG Pathway Enrichment - Upregulated Genes"
)

plot_enrichment(
    kegg_down,
    file.path(
        plots_dir,
        "KEGG_downregulated_dotplot.png"
    ),
    "KEGG Pathway Enrichment - Downregulated Genes"
)

if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {

    gsea_plot <- dotplot(
        gsea_go,
        showCategory = 15
    ) +
        ggtitle(
            "GSEA - GO Biological Process"
        )

    ggsave(
        filename = file.path(
            plots_dir,
            "GSEA_GO_BP_dotplot.png"
        ),
        plot = gsea_plot,
        width = 10,
        height = 7,
        dpi = 300
    )

    cat(
        "Created:",
        file.path(
            plots_dir,
            "GSEA_GO_BP_dotplot.png"
        ),
        "\n"
    )
}

# =========================================================
# Save summary
# =========================================================

summary_lines <- c(
    "RNAFlowX Functional Enrichment Summary",
    "=======================================",
    paste(
        "DESeq2 genes:",
        nrow(res)
    ),
    paste(
        "Mapped genes:",
        nrow(mapped_res)
    ),
    paste(
        "Significant genes:",
        nrow(significant)
    ),
    paste(
        "Upregulated genes:",
        nrow(upregulated)
    ),
    paste(
        "Downregulated genes:",
        nrow(downregulated)
    ),
    paste(
        "GSEA ranked genes:",
        length(gene_list)
    )
)

writeLines(
    summary_lines,
    file.path(
        outdir,
        "enrichment_summary.txt"
    )
)

cat("\n====================================\n")
cat("ENRICHMENT COMPLETED SUCCESSFULLY\n")
cat("====================================\n")
