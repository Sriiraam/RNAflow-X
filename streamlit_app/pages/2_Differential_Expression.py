import streamlit as st
import pandas as pd
from pathlib import Path

# ============================================================
# PAGE CONFIG
# ============================================================
st.set_page_config(
    page_title="RNAFlowX | Differential Expression",
    layout="wide"
)

# ============================================================
# PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULT_DIR = (
    PROJECT_ROOT
    / "results"
    / "differential_expression"
    / "deseq2_results"
)

DE_FILE = RESULT_DIR / "differential_expression.csv"
PCA_FILE = RESULT_DIR / "PCA.png"
VOLCANO_FILE = RESULT_DIR / "volcano_plot.png"
MA_FILE = RESULT_DIR / "MA_plot.png"

# ============================================================
# CSS
# ============================================================
st.markdown(
    """
<style>

/* ---------- HERO ---------- */
.de-hero {
    position: relative;
    overflow: hidden;
    padding: 44px 46px 48px;
    margin-bottom: 30px;

    border-radius: 26px;

    background:
        linear-gradient(
            135deg,
            #062A67 0%,
            #0B4BC4 55%,
            #1769E0 100%
        );

    box-shadow:
        0 18px 38px rgba(15, 74, 180, 0.20);

    color: white;
}

.de-hero::after {
    content: "";
    position: absolute;

    width: 420px;
    height: 420px;

    right: -150px;
    bottom: -240px;

    border-radius: 50%;

    background:
        radial-gradient(
            circle,
            rgba(109, 227, 255, 0.22),
            transparent 68%
        );
}

.de-hero h1 {
    margin: 8px 0 10px;
    font-size: 48px;
    font-weight: 800;
    letter-spacing: -0.5px;
}

.de-hero h3 {
    margin: 0 0 18px;
    font-size: 22px;
    color: #72E6FF;
    font-weight: 600;
}

.de-hero p {
    max-width: 920px;
    margin: 0;

    color: #E8F3FF;
    font-size: 16px;
    line-height: 1.8;
}

.plot-title{
    font-size:28px;
    font-weight:800;
    color:#1E293B;
    margin-bottom:14px;
}

.plot-caption{
    text-align:center;
    font-size:15px;
    font-weight:600;
    color:#64748B;
    margin-top:8px;
}

/* ---------- KPI ---------- */
.de-kpi {
    min-height: 145px;

    background: white;
    border-radius: 18px;

    padding: 20px;

    border-top: 6px solid;

    box-shadow:
        0 7px 20px rgba(15, 23, 42, 0.07);

    display: flex;
    flex-direction: column;
    justify-content: center;
    text-align: center;
}

.de-kpi-icon {
    font-size: 27px;
    margin-bottom: 6px;
}

.de-kpi-value {
    font-size: 37px;
    line-height: 1.1;
    font-weight: 800;
    color: #0F172A;
}

.de-kpi-label {
    margin-top: 8px;
    font-size: 14px;
    font-weight: 600;
    color: #64748B;
}

/* ---------- SECTION ---------- */
.de-section{
    margin:40px 0 20px;
    font-size:36px;
    font-weight:800;
    color:#0F3F7A;
}

/* ---------- ANALYSIS CARD ---------- */
.de-card {
    background: white;

    border: 1px solid #E5E7EB;
    border-radius: 18px;

    padding: 20px;

    box-shadow:
        0 6px 18px rgba(15, 23, 42, 0.05);
}

/* ---------- SUMMARY ---------- */
.de-summary {
    background: #EFF6FF;

    border-left: 5px solid #2563EB;
    border-radius: 14px;

    padding: 20px 22px;

    color: #334155;

    font-size: 15px;
    line-height: 1.8;
}

.de-summary strong {
    color: #1D4ED8;
}

</style>
""",
    unsafe_allow_html=True
)

# ============================================================
# LOAD DATA
# ============================================================
if not DE_FILE.exists():
    st.error(
        f"Differential-expression file not found:\n\n{DE_FILE}"
    )
    st.stop()

deg = pd.read_csv(DE_FILE, index_col=0)

# Convert Ensembl row names into a real column
deg.index.name = "Gene"
deg = deg.reset_index()

required_columns = {
    "Gene",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
}

missing_columns = required_columns - set(deg.columns)

if missing_columns:
    st.error(
        "Missing required DESeq2 columns: "
        + ", ".join(sorted(missing_columns))
    )
    st.stop()

# Clean values
deg["padj"] = pd.to_numeric(
    deg["padj"],
    errors="coerce"
).fillna(1.0)

deg["pvalue"] = pd.to_numeric(
    deg["pvalue"],
    errors="coerce"
).fillna(1.0)

deg["log2FoldChange"] = pd.to_numeric(
    deg["log2FoldChange"],
    errors="coerce"
)

deg["baseMean"] = pd.to_numeric(
    deg["baseMean"],
    errors="coerce"
)

# ============================================================
# HERO
# ============================================================
st.markdown(
    """
<div class="de-hero">

<h1>Differential Expression Analysis</h1>

<h3>
DESeq2 • Statistical Transcriptomics • PCA • Volcano • MA Plot
</h3>

<p>
This module evaluates transcriptional changes between control and
PFOS-treated samples using DESeq2. Statistical significance is assessed
using Benjamini-Hochberg adjusted p-values together with configurable
fold-change thresholds.
</p>

</div>
""",
    unsafe_allow_html=True
)

# ============================================================
# ANALYSIS THRESHOLDS
# ============================================================
with st.expander(
    "Analysis thresholds",
    expanded=False
):
    threshold_col1, threshold_col2 = st.columns(2)

    with threshold_col1:
        padj_threshold = st.number_input(
            "Adjusted p-value threshold",
            min_value=0.001,
            max_value=0.10,
            value=0.05,
            step=0.01,
            format="%.3f"
        )

    with threshold_col2:
        lfc_threshold = st.number_input(
            "|log2 Fold Change| threshold",
            min_value=0.0,
            max_value=5.0,
            value=1.0,
            step=0.25
        )

# ============================================================
# FILTER RESULTS
# ============================================================
significant = deg[
    (deg["padj"] < padj_threshold)
    & (deg["log2FoldChange"].abs() >= lfc_threshold)
].copy()

upregulated = significant[
    significant["log2FoldChange"] >= lfc_threshold
].copy()

downregulated = significant[
    significant["log2FoldChange"] <= -lfc_threshold
].copy()

# ============================================================
# KPI CARDS
# ============================================================
k1, k2, k3, k4 = st.columns(
    4,
    gap="medium"
)

with k1:
    st.markdown(
        f"""
<div class="de-kpi" style="border-color:#2563EB;">
<div class="de-kpi-icon">🧬</div>
<div class="de-kpi-value">{len(deg):,}</div>
<div class="de-kpi-label">Genes Tested</div>
</div>
""",
        unsafe_allow_html=True
    )

with k2:
    st.markdown(
        f"""
<div class="de-kpi" style="border-color:#16A34A;">
<div class="de-kpi-icon">◆</div>
<div class="de-kpi-value">{len(significant):,}</div>
<div class="de-kpi-label">Significant DEGs</div>
</div>
""",
        unsafe_allow_html=True
    )

with k3:
    st.markdown(
        f"""
<div class="de-kpi" style="border-color:#7C3AED;">
<div class="de-kpi-icon">▲</div>
<div class="de-kpi-value">{len(upregulated):,}</div>
<div class="de-kpi-label">Upregulated</div>
</div>
""",
        unsafe_allow_html=True
    )

with k4:
    st.markdown(
        f"""
<div class="de-kpi" style="border-color:#EA580C;">
<div class="de-kpi-icon">▼</div>
<div class="de-kpi-value">{len(downregulated):,}</div>
<div class="de-kpi-label">Downregulated</div>
</div>
""",
        unsafe_allow_html=True
    )

# ============================================================
# PCA + VOLCANO
# ============================================================

st.markdown(
    '<div class="de-section">📊 Expression Landscape</div>',
    unsafe_allow_html=True
)

plot1, plot2 = st.columns(2, gap="large")

with plot1:
    st.markdown(
        '<div class="plot-title">Principal Component Analysis</div>',
        unsafe_allow_html=True
    )

    if PCA_FILE.exists():
        st.image(
            str(PCA_FILE),
            width="stretch"
        )

        st.markdown(
            '<div class="plot-caption">PCA of normalized gene-expression profiles</div>',
            unsafe_allow_html=True
        )

with plot2:
    st.markdown(
        '<div class="plot-title">Volcano Plot</div>',
        unsafe_allow_html=True
    )

    if VOLCANO_FILE.exists():
        st.image(
            str(VOLCANO_FILE),
            width="stretch"
        )

        st.markdown(
            '<div class="plot-caption">Effect size versus statistical significance</div>',
            unsafe_allow_html=True
        )

# ============================================================
# MA PLOT
# ============================================================

st.markdown(
    '<div class="de-section">📈 Expression Effect Distribution</div>',
    unsafe_allow_html=True
)

ma_left, ma_right = st.columns([1.45, 0.85], gap="large")

with ma_left:
    st.markdown(
        '<div class="plot-title">MA Plot</div>',
        unsafe_allow_html=True
    )

    if MA_FILE.exists():
        st.image(
            str(MA_FILE),
            width="stretch"
        )

        st.markdown(
            '<div class="plot-caption">Mean normalized expression versus log2 fold change</div>',
            unsafe_allow_html=True
        )

with ma_right:
    st.markdown(
        f"""
<div class="de-summary">
<strong>Current Statistical Definition</strong><br><br>

Adjusted p-value &lt; <strong>{padj_threshold}</strong><br><br>

Absolute log2 fold change ≥ <strong>{lfc_threshold}</strong><br><br>

<strong>{len(significant):,}</strong> genes satisfy both criteria.
</div>
""",
        unsafe_allow_html=True
    )

# ============================================================
# SIGNIFICANT DEG TABLE
# ============================================================
st.markdown(
    '<div class="de-section">🧬 Significant Gene Explorer</div>',
    unsafe_allow_html=True
)

search_col, regulation_col = st.columns(
    [2, 1]
)

with search_col:
    gene_search = st.text_input(
        "Search Ensembl Gene ID",
        placeholder="Example: ENSG00000101255"
    )

with regulation_col:
    regulation_filter = st.selectbox(
        "Regulation",
        [
            "All significant",
            "Upregulated",
            "Downregulated"
        ]
    )

if regulation_filter == "Upregulated":
    table = upregulated.copy()

elif regulation_filter == "Downregulated":
    table = downregulated.copy()

else:
    table = significant.copy()

if gene_search:
    table = table[
        table["Gene"]
        .astype(str)
        .str.contains(
            gene_search,
            case=False,
            na=False
        )
    ]

table = table.sort_values(
    "padj",
    ascending=True
)

display_table = table[
    [
        "Gene",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
    ]
].copy()

display_table["baseMean"] = (
    display_table["baseMean"]
    .round(1)
)

display_table["log2FoldChange"] = (
    display_table["log2FoldChange"]
    .round(2)
)

styled_table = (
    display_table.style
    .format({
        "baseMean": "{:,.1f}",
        "log2FoldChange": "{:.2f}",
        "pvalue": "{:.2e}",
        "padj": "{:.2e}"
    })
)


st.dataframe(
    styled_table,
    width="stretch",
    hide_index=True,
    height=430
)

st.caption(
    f"Displaying {len(display_table):,} genes "
    f"after the current filters."
)

# ============================================================
# TOP DIFFERENTIALLY EXPRESSED GENES
# ============================================================
st.markdown(
    '<div class="de-section">🏆 Strongest Differential Signals</div>',
    unsafe_allow_html=True
)

top_left, top_right = st.columns(
    2,
    gap="large"
)

with top_left:
    st.markdown("### Top Upregulated")

    top_up = (
        upregulated
        .sort_values(
            "log2FoldChange",
            ascending=False
        )
        .head(10)
    )

    st.dataframe(
        top_up[
            [
                "Gene",
                "log2FoldChange",
                "padj"
            ]
        ].style.format({
            "log2FoldChange": "{:.2f}",
            "padj": "{:.2e}"
        }),
        width="stretch",
        hide_index=True
    )

with top_right:
    st.markdown("### Top Downregulated")

    top_down = (
        downregulated
        .sort_values(
            "log2FoldChange",
            ascending=True
        )
        .head(10)
    )

    st.dataframe(
        top_down[
            [
                "Gene",
                "log2FoldChange",
                "padj"
            ]
        ].style.format({
            "log2FoldChange": "{:.2f}",
            "padj": "{:.2e}"
        }),
        width="stretch",
        hide_index=True
    )

# ============================================================
# ANALYSIS SUMMARY
# ============================================================
st.markdown(
    '<div class="de-section">🔬 Statistical Interpretation</div>',
    unsafe_allow_html=True
)

st.markdown(
    f"""
<div class="de-summary">

The DESeq2 analysis tested
<strong>{len(deg):,}</strong> genes.

Using an adjusted p-value threshold of
<strong>{padj_threshold}</strong>
and an absolute log2 fold-change threshold of
<strong>{lfc_threshold}</strong>:

<br><br>

• <strong>{len(significant):,}</strong>
genes are classified as significantly differentially expressed.<br>

• <strong>{len(upregulated):,}</strong>
genes are classified as upregulated.<br>

• <strong>{len(downregulated):,}</strong>
genes are classified as downregulated.<br><br>

These significant genes provide the input for downstream
GO, KEGG and GSEA pathway-level interpretation.

</div>
""",
    unsafe_allow_html=True
)

st.write("")

st.caption(
    "Method: DESeq2 • Multiple-testing correction: "
    "Benjamini-Hochberg FDR"
)
