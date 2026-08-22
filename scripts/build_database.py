from pathlib import Path
import sqlite3
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DB_DIR = ROOT / "database"
DB_PATH = DB_DIR / "rnaflowx.db"

DB_DIR.mkdir(parents=True, exist_ok=True)

SOURCES = {
    "samples": ROOT / "data/metadata.csv",

    "qc_metrics":
        ROOT / "results/multiqc/multiqc_report_data/multiqc_fastqc.txt",

    "differential_expression":
        ROOT / "results/differential_expression/deseq2_results/differential_expression.csv",

    "go_bp_all":
        ROOT / "results/enrichment/GO/BP/all_significant.csv",

    "go_bp_up":
        ROOT / "results/enrichment/GO/BP/upregulated.csv",

    "go_bp_down":
        ROOT / "results/enrichment/GO/BP/downregulated.csv",

    "kegg_all":
        ROOT / "results/enrichment/KEGG/all_significant.csv",

    "kegg_up":
        ROOT / "results/enrichment/KEGG/upregulated.csv",

    "kegg_down":
        ROOT / "results/enrichment/KEGG/downregulated.csv",

    "gsea_go_bp":
        ROOT / "results/enrichment/GSEA/GO_BP_GSEA.csv",

    "go_mf_all":
        ROOT / "results/enrichment/GO/MF/all_significant.csv",

    "go_mf_up":
        ROOT / "results/enrichment/GO/MF/upregulated.csv",

    "go_mf_down":
        ROOT / "results/enrichment/GO/MF/downregulated.csv",

    "go_cc_all":
        ROOT / "results/enrichment/GO/CC/all_significant.csv",

    "go_cc_up":
        ROOT / "results/enrichment/GO/CC/upregulated.csv",

    "go_cc_down":
        ROOT / "results/enrichment/GO/CC/downregulated.csv",
}


def load_csv(conn, table_name, path):
    if not path.exists():
        print(f"[SKIP] {table_name}: {path} not found")
        return

    if path.suffix == ".txt":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)

    # Preserve gene identifiers in DESeq2 results.
    if table_name == "differential_expression":
        first_col = df.columns[0]

        if first_col.startswith("Unnamed"):
            df = df.rename(columns={first_col: "Gene"})
    else:
        df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

    df.to_sql(
        table_name,
        conn,
        if_exists="replace",
        index=False
    )

    print(f"[OK] {table_name}: {len(df):,} rows")

def create_indexes(conn):
    commands = [
        """
        CREATE INDEX IF NOT EXISTS idx_de_padj
        ON differential_expression(padj)
        """,

        """
        CREATE INDEX IF NOT EXISTS idx_de_log2fc
        ON differential_expression(log2FoldChange)
        """,

        """
        CREATE INDEX IF NOT EXISTS idx_go_bp_padj
        ON go_bp_all("p.adjust")
        """,

        """
        CREATE INDEX IF NOT EXISTS idx_kegg_padj
        ON kegg_all("p.adjust")
        """,

        """
        CREATE INDEX IF NOT EXISTS idx_gsea_nes
        ON gsea_go_bp(NES)
        """
    ]

    for command in commands:
        try:
            conn.execute(command)
        except sqlite3.OperationalError as exc:
            print(f"[INDEX SKIP] {exc}")


def main():
    if DB_PATH.exists():
        DB_PATH.unlink()

    print(f"Building RNAFlowX database:")
    print(DB_PATH)
    print()

    with sqlite3.connect(DB_PATH) as conn:

        for table_name, path in SOURCES.items():
            load_csv(conn, table_name, path)

        create_indexes(conn)

        conn.commit()

    print()
    print("RNAFlowX SQLite database created successfully.")


if __name__ == "__main__":
    main()
