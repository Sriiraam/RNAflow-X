process QC_SUMMARY {

    label 'process_low'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_files

    output:
    path "qc_summary.csv", emit: csv
    path "qc_summary.json", emit: json
    path "qc_status.txt", emit: status

    script:
    """
    python3 ${projectDir}/bin/summarize_qc.py
    """
}
