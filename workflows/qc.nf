workflow QC {

    take:
    samples

    main:
    FASTQC(samples)

    emit:
    fastqc = FASTQC.out

}
