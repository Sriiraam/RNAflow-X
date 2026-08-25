from pathlib import Path
import sqlite3
import pandas as pd

APP_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = APP_ROOT / "data" / "rnaflowx.db"


def get_connection():
    """Return a read-only connection to the packaged RNAFlowX database."""
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"RNAFlowX database not found: {DB_PATH}"
        )

    return sqlite3.connect(
        f"file:{DB_PATH}?mode=ro",
        uri=True
    )


def query(sql, params=None):
    with get_connection() as conn:
        return pd.read_sql_query(
            sql,
            conn,
            params=params or ()
        )


def get_samples():
    return query("SELECT * FROM samples")


def get_qc_metrics():
    return query("SELECT * FROM qc_metrics")


def get_differential_expression():
    return query("""
        SELECT *
        FROM differential_expression
        ORDER BY padj ASC
    """)


def get_significant_genes(padj=0.05, abs_log2fc=1.0):
    return query("""
        SELECT *
        FROM differential_expression
        WHERE padj IS NOT NULL
          AND padj < ?
          AND ABS(log2FoldChange) >= ?
        ORDER BY padj ASC
    """, (padj, abs_log2fc))


def get_go(go_type="BP", category="all"):
    go_type = go_type.lower()

    tables = {
        ("bp", "all"): "go_bp_all",
        ("bp", "up"): "go_bp_up",
        ("bp", "down"): "go_bp_down",
        ("mf", "all"): "go_mf_all",
        ("mf", "up"): "go_mf_up",
        ("mf", "down"): "go_mf_down",
        ("cc", "all"): "go_cc_all",
        ("cc", "up"): "go_cc_up",
        ("cc", "down"): "go_cc_down",
    }

    table = tables[(go_type, category)]

    return query(
        f'SELECT * FROM "{table}" ORDER BY "p.adjust" ASC'
    )


def get_kegg(category="all"):
    tables = {
        "all": "kegg_all",
        "up": "kegg_up",
        "down": "kegg_down"
    }

    table = tables.get(category, "kegg_all")

    return query(
        f'SELECT * FROM "{table}" ORDER BY "p.adjust" ASC'
    )


def get_gsea():
    return query("""
        SELECT *
        FROM gsea_go_bp
        ORDER BY ABS(NES) DESC
    """)
