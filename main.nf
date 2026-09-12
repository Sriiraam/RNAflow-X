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


    /*
     * Automatic run provenance
     */
    workflow.onComplete = {

        def provenanceDir = file("${params.outdir}/pipeline_info")
        provenanceDir.mkdirs()

        /*
         * Capture the final resolved parameter values used by this run.
         * Convert values to strings so Paths and Nextflow-specific objects
         * remain JSON serializable.
         */
        def resolvedParams = params.collectEntries { key, value ->
            [(key.toString()): value == null ? null : value.toString()]
        }

        def paramsJson = groovy.json.JsonOutput.toJson(resolvedParams)

        def cmd = [
            'python3',
            "${projectDir}/bin/write_provenance.py",
            '--output', "${params.outdir}/pipeline_info/provenance.json",
            '--pipeline-version', params.version.toString(),
            '--params-json', paramsJson,
            '--run-name', workflow.runName.toString(),
            '--session-id', workflow.sessionId.toString(),
            '--start', workflow.start.toString(),
            '--complete', workflow.complete.toString(),
            '--duration', workflow.duration.toString(),
            '--success', workflow.success.toString(),
            '--profile', workflow.profile ?: '',
            '--command-line', workflow.commandLine ?: '',
            '--nextflow-version', nextflow.version.toString(),
            '--container-image', params.container_image.toString(),
            '--samplesheet', params.samplesheet.toString(),
            '--salmon-index', params.salmon_index.toString(),
            '--tx2gene', params.tx2gene.toString(),
            '--metadata', params.metadata.toString()
        ]

        /*
         * ProcessBuilder requires every command argument
         * to be a Java String.
         */
        def safeCmd = cmd.collect { value ->
            if (value == null) {
                return ''
            }

            if (value instanceof Collection) {
                return value.join(',')
            }

            return value.toString()
        }

        def process = new ProcessBuilder(safeCmd)
            .directory(new File(projectDir.toString()))
            .inheritIO()
            .start()

        process.waitFor()
    }
}
