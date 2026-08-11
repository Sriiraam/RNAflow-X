process DESEQ2 {

    tag "Differential Expression"

    label 'process_medium'

    publishDir "${params.outdir}/differential_expression", mode: 'copy'

    input:
    path count_matrix
    path metadata

    output:
    path "deseq2_results", emit: results
    path "versions.yml", emit: versions

    script:
    """
    export R_LIBS_USER="\$HOME/R/library"

    mkdir -p deseq2_results

    Rscript ${projectDir}/bin/run_deseq2.R \
        --counts ${count_matrix} \
        --metadata ${metadata} \
        --outdir deseq2_results

    Rscript -e 'cat(
        "DESeq2:\\n    version: ",
        as.character(packageVersion("DESeq2")),
        "\\n",
        file = "versions.yml"
    )'
    """
}