process ENRICHMENT {

    tag "Functional enrichment"

    label 'process_medium'

    publishDir "${params.outdir}/functional_enrichment", mode: 'copy'

    input:
    path deseq2_results

    output:
    path "enrichment_results", emit: results
    path "versions.yml", emit: versions

    script:
    """
    export R_LIBS_USER="\$HOME/R/library"

    mkdir -p enrichment_results

    Rscript ${projectDir}/bin/run_enrichment.R \
        ${deseq2_results}/differential_expression.csv \
        enrichment_results

    Rscript -e 'cat(
        "clusterProfiler:\\n    version: ",
        as.character(packageVersion("clusterProfiler")),
        "\\norg.Hs.eg.db:\\n    version: ",
        as.character(packageVersion("org.Hs.eg.db")),
        "\\nenrichplot:\\n    version: ",
        as.character(packageVersion("enrichplot")),
        "\\nggplot2:\\n    version: ",
        as.character(packageVersion("ggplot2")),
        "\\n",
        file = "versions.yml"
    )'
    """
}
