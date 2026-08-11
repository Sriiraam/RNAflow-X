process TXIMPORT {

    tag "Gene-level count matrix"

    label 'process_medium'

    publishDir "${params.outdir}/counting", mode: 'copy'

    input:
    path quant_dirs
    path tx2gene
    path metadata

    output:
    path "count_matrix.csv", emit: counts
    path "tximport_summary.txt", emit: summary

    script:
    """
    export R_LIBS_USER="\$HOME/R/library"

    Rscript ${projectDir}/bin/run_tximport.R \
        --tx2gene ${tx2gene} \
        --metadata ${metadata} \
        --counts count_matrix.csv \
        --summary tximport_summary.txt
    """
}