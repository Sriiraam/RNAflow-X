from pathlib import Path
import pandas as pd

RESULTS = Path("tests/results_ci")

de_file = (
    RESULTS
    / "differential_expression"
    / "deseq2_results"
    / "differential_expression.csv"
)

qc_file = RESULTS / "qc" / "qc_summary.csv"


def main():
    assert de_file.exists(), f"Missing: {de_file}"
    assert qc_file.exists(), f"Missing: {qc_file}"

    de = pd.read_csv(de_file)

    assert "gene_id" in de.columns, (
        "DESeq2 output must contain an explicit gene_id column"
    )

    de = de.set_index("gene_id")

    assert "GENE001.1" in de.index
    assert "GENE006.1" in de.index

    gene1_fc = de.loc["GENE001.1", "log2FoldChange"]
    gene6_fc = de.loc["GENE006.1", "log2FoldChange"]

    assert gene1_fc > 1, (
        f"GENE001.1 expected positive log2FC, got {gene1_fc}"
    )

    assert gene6_fc < -1, (
        f"GENE006.1 expected negative log2FC, got {gene6_fc}"
    )

    qc = pd.read_csv(qc_file)

    assert len(qc) == 4, (
        f"Expected 4 QC samples, found {len(qc)}"
    )

    failed = qc[
        qc["overall_status"] == "FAIL"
    ]

    assert failed.empty, (
        "CI QC contains FAIL samples:\n"
        + failed.to_string(index=False)
    )

    expected_samples = {
        "control_rep1",
        "control_rep2",
        "treated_rep1",
        "treated_rep2",
    }

    observed_samples = set(qc["sample_id"])

    assert observed_samples == expected_samples, (
        f"Unexpected samples: {observed_samples}"
    )

    print("✅ RNAFlowX miniature end-to-end CI passed")
    print(f"GENE001.1 log2FC: {gene1_fc:.3f}")
    print(f"GENE006.1 log2FC: {gene6_fc:.3f}")
    print("QC samples: 4")
    print("QC failures: 0")


if __name__ == "__main__":
    main()
