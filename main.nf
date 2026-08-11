nextflow.enable.dsl=2

include {
    RNAFLOWX
} from './workflows/rnaseq'


workflow {

    /*
     * RNAFlowX
     * Bulk RNA-seq analysis pipeline
     */

    Channel
        .fromPath(
            params.samplesheet,
            checkIfExists: true
        )
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true)
            )
        }
        .set { sample_reads }


    /*
     * Salmon transcriptome index
     */
    salmon_index = file(
        params.salmon_index,
        checkIfExists: true
    )


    /*
     * tx2gene mapping
     */
    tx2gene = file(
        params.tx2gene,
        checkIfExists: true
    )


    /*
     * Sample metadata
     */
    metadata = file(
        params.metadata,
        checkIfExists: true
    )


    /*
     * Complete RNAFlowX workflow
     */
    RNAFLOWX(
        sample_reads,
        salmon_index,
        tx2gene,
        metadata
    )
}