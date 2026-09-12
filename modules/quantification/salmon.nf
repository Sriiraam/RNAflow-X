process SALMON_QUANT {

    tag "$sample_id"

    label 'process_low'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path salmon_index

    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant
    tuple val(sample_id), path("${sample_id}/aux_info/meta_info.json"), emit: meta
    path "${sample_id}.log", emit: logs

    script:
    """
    salmon quant \
        --index ${salmon_index} \
        --libType A \
        -1 ${read1} \
        -2 ${read2} \
        --output ${sample_id} \
        --threads ${task.cpus} \
        > ${sample_id}.log 2>&1
    """
}
