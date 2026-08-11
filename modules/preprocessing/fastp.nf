process FASTP {

    tag "$sample_id"

    label 'process_medium'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp_report.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastp_report.json"), emit: json

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_R1.fastq.gz \
        --out2 ${sample_id}_R2.fastq.gz \
        --thread ${task.cpus} \
        --html ${sample_id}_fastp_report.html \
        --json ${sample_id}_fastp_report.json \
        --detect_adapter_for_pe
    """
}