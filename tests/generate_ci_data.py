#!/usr/bin/env python3

import argparse
import gzip
import random
from pathlib import Path


DNA = "ACGT"


def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def make_transcript(rng, length=500):
    return "".join(rng.choice(DNA) for _ in range(length))


def write_fastq_record(handle, name, sequence):
    quality = "I" * len(sequence)

    handle.write(f"@{name}\n")
    handle.write(f"{sequence}\n")
    handle.write("+\n")
    handle.write(f"{quality}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir",
        default="tests/data/ci"
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    fastq_dir = outdir / "fastq"
    reference_dir = outdir / "reference"

    fastq_dir.mkdir(parents=True, exist_ok=True)
    reference_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(42)

    # -----------------------------------------------------
    # Synthetic transcriptome
    # -----------------------------------------------------

    number_of_genes = 40
    transcript_length = 500

    transcripts = {}

    for i in range(1, number_of_genes + 1):
        transcript_id = f"TX{i:03d}.1"
        gene_id = f"GENE{i:03d}.1"

        transcripts[transcript_id] = {
            "gene": gene_id,
            "sequence": make_transcript(
                rng,
                transcript_length
            )
        }

    fasta_path = reference_dir / "transcripts.fa"

    with fasta_path.open("w") as fh:
        for transcript_id, info in transcripts.items():
            fh.write(f">{transcript_id}\n")
            fh.write(info["sequence"] + "\n")

    tx2gene_path = reference_dir / "tx2gene.tsv"

    with tx2gene_path.open("w") as fh:
        for transcript_id, info in transcripts.items():
            fh.write(
                f"{transcript_id}\t{info['gene']}\n"
            )

    # -----------------------------------------------------
    # Experimental design
    # -----------------------------------------------------

    samples = {
        "control_rep1": "control",
        "control_rep2": "control",
        "treated_rep1": "treated",
        "treated_rep2": "treated",
    }

    metadata_path = outdir / "metadata.csv"

    with metadata_path.open("w") as fh:
        fh.write("sample_id,condition\n")

        for sample_id, condition in samples.items():
            fh.write(
                f"{sample_id},{condition}\n"
            )

    # -----------------------------------------------------
    # Synthetic expression model
    #
    # Genes 1-5  : strongly upregulated in treatment
    # Genes 6-10 : strongly downregulated in treatment
    # Genes 11-40: approximately stable
    # -----------------------------------------------------

    expression = {}

    for sample_id, condition in samples.items():

        sample_counts = {}

        for gene_number in range(
            1,
            number_of_genes + 1
        ):

            if 1 <= gene_number <= 5:

                if condition == "control":
                    base = 30
                else:
                    base = 180

            elif 6 <= gene_number <= 10:

                if condition == "control":
                    base = 180
                else:
                    base = 30

            else:
                base = 80

            # Deterministic biological-style replicate variation.
            #
            # The CI data must not contain nearly identical replicates,
            # because DESeq2 needs realistic gene-wise dispersion.
            expression_rng = random.Random(
                100000
                + gene_number * 100
                + list(samples).index(sample_id)
            )

            biological_jitter = expression_rng.lognormvariate(
                0.0,
                0.30
            )

            sample_scale = {
                "control_rep1": 0.92,
                "control_rep2": 1.08,
                "treated_rep1": 0.95,
                "treated_rep2": 1.05,
            }[sample_id]

            count = round(
                base
                * biological_jitter
                * sample_scale
            )

            sample_counts[gene_number] = max(
                count,
                15
            )

        expression[sample_id] = sample_counts

    # -----------------------------------------------------
    # Generate paired-end FASTQ
    # -----------------------------------------------------

    read_length = 75
    fragment_length = 200

    for sample_id in samples:

        r1_path = fastq_dir / f"{sample_id}_R1.fastq.gz"
        r2_path = fastq_dir / f"{sample_id}_R2.fastq.gz"

        read_number = 0

        with gzip.open(
            r1_path,
            "wt"
        ) as r1_fh, gzip.open(
            r2_path,
            "wt"
        ) as r2_fh:

            for gene_number in range(
                1,
                number_of_genes + 1
            ):

                transcript_id = (
                    f"TX{gene_number:03d}.1"
                )

                sequence = transcripts[
                    transcript_id
                ]["sequence"]

                n_pairs = expression[
                    sample_id
                ][gene_number]

                max_start = (
                    len(sequence)
                    - fragment_length
                )

                for pair_index in range(n_pairs):

                    read_number += 1

                    # Deterministic but distributed fragments.
                    start = (
                        pair_index * 17
                        + gene_number * 11
                    ) % (max_start + 1)

                    fragment = sequence[
                        start:
                        start + fragment_length
                    ]

                    read1 = fragment[
                        :read_length
                    ]

                    read2 = reverse_complement(
                        fragment[
                            -read_length:
                        ]
                    )

                    read_name = (
                        f"{sample_id}_"
                        f"G{gene_number:03d}_"
                        f"{read_number:06d}"
                    )

                    write_fastq_record(
                        r1_fh,
                        read_name + "/1",
                        read1
                    )

                    write_fastq_record(
                        r2_fh,
                        read_name + "/2",
                        read2
                    )

    # -----------------------------------------------------
    # Samplesheet
    # -----------------------------------------------------

    samplesheet_path = outdir / "samplesheet.csv"

    with samplesheet_path.open("w") as fh:

        fh.write(
            "sample_id,condition,"
            "fastq_1,fastq_2\n"
        )

        for sample_id, condition in samples.items():

            fh.write(
                f"{sample_id},"
                f"{condition},"
                f"{fastq_dir}/{sample_id}_R1.fastq.gz,"
                f"{fastq_dir}/{sample_id}_R2.fastq.gz\n"
            )

    print("Synthetic RNAFlowX CI dataset generated")
    print(f"Output: {outdir}")
    print(f"Samples: {len(samples)}")
    print(f"Genes: {number_of_genes}")


if __name__ == "__main__":
    main()
