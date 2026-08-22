import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(
    page_title="RNAFlowX | GSEA",
    layout="wide"
)

# ============================================================
# PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]

GSEA_FILE = (
    PROJECT_ROOT
    / "results"
    / "enrichment"
    / "GSEA"
    / "GO_BP_GSEA.csv"
)

GSEA_PLOT = (
    PROJECT_ROOT
    / "results"
    / "enrichment"
    / "plots"
    / "GSEA_GO_BP_dotplot.png"
)

# ============================================================
# CSS
# ============================================================
st.markdown("""
<style>

.gsea-hero{
    padding:46px;
    border-radius:26px;
    background:linear-gradient(135deg,#062A67 0%,#0E56D8 100%);
    color:white;
    box-shadow:0 18px 36px rgba(15,74,180,.20);
    margin-bottom:30px;
}

.gsea-hero h1{
    margin:0 0 12px;
    font-size:48px;
    font-weight:800;
}

.gsea-hero h3{
    margin:0 0 18px;
    color:#75E6FF;
    font-size:22px;
}

.gsea-hero p{
    max-width:900px;
    color:#E8F3FF;
    font-size:16px;
    line-height:1.8;
}

.kpi{
    background:white;
    border-radius:18px;
    padding:20px;
    text-align:center;
    border-top:6px solid;
    box-shadow:0 7px 20px rgba(15,23,42,.07);
    min-height:140px;
}

.kpi-value{
    font-size:38px;
    font-weight:800;
    color:#0F172A;
}

.kpi-label{
    margin-top:8px;
    font-size:14px;
    font-weight:600;
    color:#64748B;
}

.section-title{
    margin:34px 0 16px;
    font-size:30px;
    font-weight:800;
    color:#153F73;
}

.info-box{
    background:#EFF6FF;
    border-left:5px solid #2563EB;
    padding:18px;
    border-radius:14px;
    color:#334155;
    line-height:1.8;
}

</style>
""", unsafe_allow_html=True)

# ============================================================
# LOAD DATA
# ============================================================
if not GSEA_FILE.exists():
    st.error(f"GSEA file not found: {GSEA_FILE}")
    st.stop()

gsea = pd.read_csv(GSEA_FILE)

required = {
    "ID",
    "Description",
    "setSize",
    "enrichmentScore",
    "NES",
    "pvalue",
    "p.adjust",
    "qvalue",
    "rank"
}

missing = required - set(gsea.columns)

if missing:
    st.error(
        "Missing required GSEA columns: "
        + ", ".join(sorted(missing))
    )
    st.stop()

gsea["NES"] = pd.to_numeric(gsea["NES"], errors="coerce")
gsea["p.adjust"] = pd.to_numeric(gsea["p.adjust"], errors="coerce")
gsea["qvalue"] = pd.to_numeric(gsea["qvalue"], errors="coerce")

significant = gsea[gsea["p.adjust"] < 0.05].copy()

positive = significant[
    significant["NES"] > 0
].copy()

negative = significant[
    significant["NES"] < 0
].copy()

# ============================================================
# HERO
# ============================================================
st.markdown("""
<div class="gsea-hero">

<h1>Gene Set Enrichment Analysis</h1>

<h3>
GSEA • Normalized Enrichment Score • GO Biological Process
</h3>

<p>
Gene Set Enrichment Analysis evaluates coordinated expression shifts across
biological pathways without requiring an arbitrary differential-expression cutoff.
RNAFlowX uses ranked gene-level statistics to identify positively and negatively
enriched Gene Ontology Biological Process gene sets.
</p>

</div>
""", unsafe_allow_html=True)

# ============================================================
# KPI
# ============================================================
k1, k2, k3, k4 = st.columns(4)

with k1:
    st.markdown(f"""
<div class="kpi" style="border-color:#2563EB">
<div class="kpi-value">{len(gsea):,}</div>
<div class="kpi-label">Gene Sets Tested</div>
</div>
""", unsafe_allow_html=True)

with k2:
    st.markdown(f"""
<div class="kpi" style="border-color:#16A34A">
<div class="kpi-value">{len(significant):,}</div>
<div class="kpi-label">Significant Gene Sets</div>
</div>
""", unsafe_allow_html=True)

with k3:
    st.markdown(f"""
<div class="kpi" style="border-color:#7C3AED">
<div class="kpi-value">{len(positive):,}</div>
<div class="kpi-label">Positive NES</div>
</div>
""", unsafe_allow_html=True)

with k4:
    st.markdown(f"""
<div class="kpi" style="border-color:#EA580C">
<div class="kpi-value">{len(negative):,}</div>
<div class="kpi-label">Negative NES</div>
</div>
""", unsafe_allow_html=True)

# ============================================================
# GSEA VISUALIZATION
# ============================================================
st.markdown(
    '<div class="section-title">📊 GSEA Enrichment Landscape</div>',
    unsafe_allow_html=True
)

left, right = st.columns([1.35, 0.9], gap="large")

with left:
    if GSEA_PLOT.exists():
        st.image(
            str(GSEA_PLOT),
            caption="GO Biological Process GSEA dotplot",
            width="stretch"
        )
    else:
        st.warning("GSEA dotplot not found.")

with right:
    strongest_positive = (
        positive.sort_values("NES", ascending=False).iloc[0]
        if not positive.empty
        else None
    )

    strongest_negative = (
        negative.sort_values("NES", ascending=True).iloc[0]
        if not negative.empty
        else None
    )

    st.markdown("### Enrichment Summary")

    if strongest_positive is not None:
        st.success(
            f"Strongest positive enrichment:\n\n"
            f"**{strongest_positive['Description']}**\n\n"
            f"NES = {strongest_positive['NES']:.2f}"
        )

    if strongest_negative is not None:
        st.error(
            f"Strongest negative enrichment:\n\n"
            f"**{strongest_negative['Description']}**\n\n"
            f"NES = {strongest_negative['NES']:.2f}"
        )

# ============================================================
# GSEA EXPLORER
# ============================================================
st.markdown(
    '<div class="section-title">🧬 GSEA Gene-Set Explorer</div>',
    unsafe_allow_html=True
)

search_col, direction_col = st.columns([2, 1])

with search_col:
    search = st.text_input(
        "Search biological process",
        placeholder="Example: cholesterol biosynthetic process"
    )

with direction_col:
    direction = st.selectbox(
        "Enrichment direction",
        [
            "All Significant",
            "Positive NES",
            "Negative NES"
        ]
    )

if direction == "Positive NES":
    table = positive.copy()

elif direction == "Negative NES":
    table = negative.copy()

else:
    table = significant.copy()

if search:
    table = table[
        table["Description"]
        .astype(str)
        .str.contains(
            search,
            case=False,
            na=False
        )
    ]

table = table.sort_values(
    "p.adjust",
    ascending=True
)

display_table = table[
    [
        "ID",
        "Description",
        "setSize",
        "NES",
        "enrichmentScore",
        "pvalue",
        "p.adjust",
        "qvalue",
        "rank"
    ]
].copy()

styled_table = (
    display_table.style
    .format({
        "NES": "{:.2f}",
        "enrichmentScore": "{:.3f}",
        "pvalue": "{:.2e}",
        "p.adjust": "{:.2e}",
        "qvalue": "{:.2e}"
    })
)

st.dataframe(
    styled_table,
    width="stretch",
    hide_index=True,
    height=430
)

st.caption(
    f"{len(display_table):,} significant gene sets displayed."
)

# ============================================================
# TOP SIGNALS
# ============================================================
st.markdown(
    '<div class="section-title">🏆 Strongest Enrichment Signals</div>',
    unsafe_allow_html=True
)

top_left, top_right = st.columns(2, gap="large")

with top_left:
    st.markdown("### Top Positive NES")

    top_positive = (
        positive
        .sort_values("NES", ascending=False)
        .head(10)
    )

    if not top_positive.empty:
        st.dataframe(
            top_positive[
                [
                    "Description",
                    "setSize",
                    "NES",
                    "p.adjust"
                ]
            ].style.format({
                "NES": "{:.2f}",
                "p.adjust": "{:.2e}"
            }),
            width="stretch",
            hide_index=True
        )

with top_right:
    st.markdown("### Top Negative NES")

    top_negative = (
        negative
        .sort_values("NES", ascending=True)
        .head(10)
    )

    if not top_negative.empty:
        st.dataframe(
            top_negative[
                [
                    "Description",
                    "setSize",
                    "NES",
                    "p.adjust"
                ]
            ].style.format({
                "NES": "{:.2f}",
                "p.adjust": "{:.2e}"
            }),
            width="stretch",
            hide_index=True
        )

# ============================================================
# CORE ENRICHMENT
# ============================================================
st.markdown(
    '<div class="section-title">🔎 Leading-Edge Signal</div>',
    unsafe_allow_html=True
)

leading_table = significant[
    [
        "Description",
        "NES",
        "leading_edge",
        "core_enrichment"
    ]
].copy()

leading_table = leading_table.sort_values(
    "NES",
    key=lambda s: s.abs(),
    ascending=False
).head(20)

st.dataframe(
    leading_table.style.format({
        "NES": "{:.2f}"
    }),
    width="stretch",
    hide_index=True,
    height=360
)

# ============================================================
# INTERPRETATION
# ============================================================
st.markdown(
    '<div class="section-title">🔬 Biological Interpretation</div>',
    unsafe_allow_html=True
)

st.markdown("""
<div class="info-box">

<strong>How to interpret the GSEA results</strong><br><br>

Positive NES values indicate gene sets whose members tend to occur toward the
upregulated end of the ranked gene list.

Negative NES values indicate gene sets concentrated toward the downregulated end.

In the current RNAFlowX results, strongly negative enrichment is observed in
sterol and cholesterol biosynthetic processes, while positive enrichment includes
coordinated pathways such as cytoplasmic translation.

GSEA therefore complements conventional over-representation analysis by detecting
coordinated biological changes even when individual genes do not cross a strict
differential-expression cutoff.

</div>
""", unsafe_allow_html=True)

st.write("")

st.caption(
    "GSEA source: clusterProfiler • GO Biological Process • "
    "Significance criterion: adjusted p-value < 0.05"
)
