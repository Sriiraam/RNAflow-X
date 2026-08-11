nextflow.enable.dsl=2

include {
    SALMON_QUANT
} from '../modules/quantification/salmon'

workflow QUANTIFICATION {

    take:
    reads
    salmon_index

    main:
    SALMON_QUANT(
        reads,
        salmon_index
    )

    emit:
    quant = SALMON_QUANT.out.quant
    logs  = SALMON_QUANT.out.logs
}
