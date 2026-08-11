process MULTIQC {

    tag "RNAFlowX MultiQC"

    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

   output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    path "versions.yml", emit: versions

    script:
    """
    multiqc . \
        --force \
        --filename multiqc_report.html \
        --outdir .

    cat > versions.yml <<EOF
MULTIQC:
    version: \$(multiqc --version | awk '{print \$3}')
EOF
    """
}