process QC_EVALUATION {

    tag "$sample_id"

    label 'process_low'

    publishDir "${params.outdir}/qc/per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(fastp_json), path(salmon_meta)

    output:
    tuple val(sample_id), path("${sample_id}_qc.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}_qc.json"), emit: json
    tuple val(sample_id), path("${sample_id}_qc_status.txt"), emit: status

    script:
    """
    python3 ${projectDir}/bin/evaluate_qc.py \
        --sample ${sample_id} \
        --fastp ${fastp_json} \
        --salmon ${salmon_meta} \
        --q30-pass ${params.qc_q30_pass} \
        --q30-warn ${params.qc_q30_warn} \
        --retention-pass ${params.qc_retention_pass} \
        --retention-warn ${params.qc_retention_warn} \
        --mapping-pass ${params.qc_mapping_pass} \
        --mapping-warn ${params.qc_mapping_warn}
    """
}
