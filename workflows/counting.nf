nextflow.enable.dsl=2

include {
    TXIMPORT
} from '../modules/counting/tximport'

workflow COUNTING {

    take:
    quant_dirs
    tx2gene
    metadata

    main:

    TXIMPORT(
        quant_dirs,
        tx2gene,
        metadata
    )

    emit:
    counts = TXIMPORT.out.counts
    summary = TXIMPORT.out.summary
}