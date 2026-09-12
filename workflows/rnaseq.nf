nextflow.enable.dsl=2

include {
    FASTQC as FASTQC_RAW
} from '../modules/qc/fastqc'

include {
    FASTQC as FASTQC_TRIMMED
} from '../modules/qc/fastqc'

include {
    FASTP
} from '../modules/preprocessing/fastp'

include {
    SALMON_QUANT
} from '../modules/quantification/salmon'

include {
    TXIMPORT
} from '../modules/counting/tximport'

include {
    DESEQ2
} from '../modules/differential_expression/deseq2'

include {
    ENRICHMENT
} from '../modules/functional_enrichment/enrichment'

include {
    MULTIQC
} from '../modules/reporting/multiqc'

include {
    QC_EVALUATION
} from '../modules/qc/qc_evaluation'

include {
    QC_SUMMARY
} from '../modules/qc/qc_summary'


workflow RNAFLOWX {

    take:
    samples
    salmon_index
    tx2gene
    metadata

    main:

    /*
     * 1. Raw-read quality control
     */
    FASTQC_RAW(samples)

    /*
     * 2. Adapter trimming and quality filtering
     */
    FASTP(samples)

    /*
     * 3. Post-trimming quality control
     */
    FASTQC_TRIMMED(FASTP.out.reads)

    /*
     * 4. Salmon transcript quantification
     */
    SALMON_QUANT(
        FASTP.out.reads,
        salmon_index
    )

    /*
     * Formal QC evaluation
     *
     * Join fastp and Salmon metrics by sample_id.
     */
    qc_inputs = FASTP.out.json
        .join(SALMON_QUANT.out.meta)

    QC_EVALUATION(qc_inputs)

    qc_tsv_files = QC_EVALUATION.out.tsv
        .map { sample_id, qc_file -> qc_file }
        .collect()

    QC_SUMMARY(qc_tsv_files)

    /*
     * 5. Collect Salmon quantification directories
     */
    salmon_quant_dirs = SALMON_QUANT.out.quant
        .map { sample_id, quant_dir -> quant_dir }
        .collect()

    /*
     * 6. Transcript-level to gene-level counts
     */
    TXIMPORT(
        salmon_quant_dirs,
        tx2gene,
        metadata
    )

    /*
     * 7. Differential expression analysis
     */
    DESEQ2(
        TXIMPORT.out.txi,
        metadata
    )

    /*
     * 8. Functional enrichment
     *    GO / KEGG / GSEA
     */
    enrichment_input = params.skip_enrichment
        ? Channel.empty()
        : DESEQ2.out.results

    ENRICHMENT(
        enrichment_input
    )

    /*
     * 9. Collect QC/reporting files for MultiQC
     *
     * FastQC:
     *   - raw FastQC zip
     *   - trimmed FastQC zip
     *
     * fastp:
     *   - HTML report
     *   - JSON report
     *
     * Salmon:
     *   - Salmon logs
     */
    multiqc_inputs = Channel
    .empty()
    .mix(FASTQC_RAW.out.zip.map { sample_id, file -> file })
    .mix(FASTQC_TRIMMED.out.zip.map { sample_id, file -> file })
    .mix(FASTP.out.html.map { sample_id, file -> file })
    .mix(FASTP.out.json.map { sample_id, file -> file })
    .mix(SALMON_QUANT.out.logs)
    .collect()

    /*
     * 10. Generate unified MultiQC report
     */
    MULTIQC(multiqc_inputs)


    emit:

    raw_qc_zip = FASTQC_RAW.out.zip
    raw_qc_html = FASTQC_RAW.out.html

    trimmed_reads = FASTP.out.reads

    trimmed_qc_zip = FASTQC_TRIMMED.out.zip
    trimmed_qc_html = FASTQC_TRIMMED.out.html

    quantification = SALMON_QUANT.out.quant
    quantification_logs = SALMON_QUANT.out.logs

    count_matrix = TXIMPORT.out.counts
    tximport_object = TXIMPORT.out.txi
    tximport_summary = TXIMPORT.out.summary

    deseq2_results = DESEQ2.out.results
    deseq2_versions = DESEQ2.out.versions

    enrichment_results = ENRICHMENT.out.results
    enrichment_versions = ENRICHMENT.out.versions

    qc_per_sample = QC_EVALUATION.out.tsv
    qc_summary_csv = QC_SUMMARY.out.csv
    qc_summary_json = QC_SUMMARY.out.json
    qc_status = QC_SUMMARY.out.status

    multiqc_report = MULTIQC.out.report
    multiqc_data = MULTIQC.out.data
    multiqc_versions = MULTIQC.out.versions
}