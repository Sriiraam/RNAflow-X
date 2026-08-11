process FASTQC {

    tag "$sample_id"

    label 'process_low'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*.zip"), emit: zip
    tuple val(sample_id), path("*.html"), emit: html

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        ${read1} ${read2}
    """
}
