from utils.database import get_go, get_kegg
import streamlit as st
import pandas as pd
from pathlib import Path

# ============================================================
# PAGE CONFIG
# ============================================================
st.set_page_config(
    page_title="RNAFlowX | Functional Enrichment",
    layout="wide"
)

# ============================================================
# PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]

ENRICH_DIR = PROJECT_ROOT / "results" / "enrichment"
GO_DIR = ENRICH_DIR / "GO"
KEGG_DIR = ENRICH_DIR / "KEGG"
PLOT_DIR = ENRICH_DIR / "plots"

# ============================================================
# CSS
# ============================================================
st.markdown("""
<style>

.enrich-hero{
    padding:46px;
    border-radius:26px;
    background:linear-gradient(135deg,#062A67 0%,#0E56D8 100%);
    color:white;
    box-shadow:0 18px 36px rgba(15,74,180,.20);
    margin-bottom:30px;
}

.enrich-hero h1{
    margin:0 0 12px;
    font-size:48px;
    font-weight:800;
}

.enrich-hero h3{
    margin:0 0 18px;
    color:#75E6FF;
    font-size:22px;
}

.enrich-hero p{
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
# HELPERS
# ============================================================
def load_csv(path):
    if not path.exists():
        return pd.DataFrame()

    return pd.read_csv(path)


def significant_count(df):
    if df.empty or "p.adjust" not in df.columns:
        return 0

    return int((df["p.adjust"] < 0.05).sum())


def prepare_table(df):
    if df.empty:
        return df

    cols = [
        c for c in [
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "Count",
            "pvalue",
            "p.adjust",
            "qvalue"
        ]
        if c in df.columns
    ]

    return df[cols].copy()


# ============================================================
# LOAD DEFAULT DATA
# ===========================================================
go_bp_all = get_go("BP", "all")
go_bp_up = get_go("BP", "up")
go_bp_down = get_go("BP", "down")

kegg_all = get_kegg("all")
kegg_up = get_kegg("up")
kegg_down = get_kegg("down")

# ============================================================
# HERO
# ============================================================
st.markdown("""
<div class="enrich-hero">

<h1>Functional Enrichment Analysis</h1>

<h3>
Gene Ontology • KEGG • Biological Interpretation
</h3>

<p>
Functional enrichment converts differentially expressed genes into
biologically interpretable pathways and processes. RNAFlowX integrates
Gene Ontology and KEGG enrichment to identify molecular functions,
cellular components, biological processes and signaling pathways
associated with PFOS exposure.
</p>

</div>
""", unsafe_allow_html=True)

# ============================================================
# KPI CARDS
# ============================================================
k1, k2, k3, k4 = st.columns(4)

with k1:
    st.markdown(f"""
<div class="kpi" style="border-color:#2563EB">
<div class="kpi-value">{significant_count(go_bp_all)}</div>
<div class="kpi-label">GO BP Terms</div>
</div>
""", unsafe_allow_html=True)

with k2:
    st.markdown(f"""
<div class="kpi" style="border-color:#16A34A">
<div class="kpi-value">{significant_count(kegg_all)}</div>
<div class="kpi-label">KEGG Pathways</div>
</div>
""", unsafe_allow_html=True)

with k3:
    st.markdown(f"""
<div class="kpi" style="border-color:#7C3AED">
<div class="kpi-value">{significant_count(go_bp_up)}</div>
<div class="kpi-label">Upregulated GO Terms</div>
</div>
""", unsafe_allow_html=True)

with k4:
    st.markdown(f"""
<div class="kpi" style="border-color:#EA580C">
<div class="kpi-value">{significant_count(go_bp_down)}</div>
<div class="kpi-label">Downregulated GO Terms</div>
</div>
""", unsafe_allow_html=True)

# ============================================================
# MAIN TABS
# ============================================================
st.markdown(
    '<div class="section-title">🧬 Enrichment Explorer</div>',
    unsafe_allow_html=True
)

tab_go, tab_kegg = st.tabs(
    ["Gene Ontology", "KEGG Pathways"]
)

# ============================================================
# GO TAB
# ============================================================
with tab_go:

    go_type = st.selectbox(
        "Gene Ontology category",
        [
            "Biological Process (BP)",
            "Molecular Function (MF)",
            "Cellular Component (CC)"
        ]
    )

    regulation = st.radio(
        "Gene set",
        [
            "All Significant",
            "Upregulated",
            "Downregulated"
        ],
        horizontal=True
    )

    go_map = {
        "Biological Process (BP)": "BP",
        "Molecular Function (MF)": "MF",
        "Cellular Component (CC)": "CC"
    }

    regulation_map = {
        "All Significant": "all",
        "Upregulated": "up",
        "Downregulated": "down"
    }

    go_code = go_map[go_type]

    go_df = get_go(
        go_code,
        regulation_map[regulation]
    )


    if go_df.empty:
        st.warning("No enrichment results available for this selection.")

    else:
        st.markdown(
            f"### {go_type} — {regulation}"
        )

        search = st.text_input(
            "Search GO term",
            placeholder="Example: vasculature development",
            key="go_search"
        )

        go_table = go_df.copy()

        if search:
            go_table = go_table[
                go_table["Description"]
                .astype(str)
                .str.contains(
                    search,
                    case=False,
                    na=False
                )
            ]

        go_table = go_table.sort_values(
            "p.adjust",
            ascending=True
        )

        display_go = prepare_table(go_table)

        st.dataframe(
            display_go.style.format({
                "pvalue": "{:.2e}",
                "p.adjust": "{:.2e}",
                "qvalue": "{:.2e}"
            }),
            width="stretch",
            hide_index=True,
            height=430
        )

        st.caption(
            f"{len(go_table)} enrichment terms displayed."
        )

# ============================================================
# KEGG TAB
# ============================================================
with tab_kegg:

    kegg_regulation = st.radio(
        "Gene set",
        [
            "All Significant",
            "Upregulated",
            "Downregulated"
        ],
        horizontal=True,
        key="kegg_regulation"
    )

    kegg_map = {
        "All Significant": "all",
        "Upregulated": "up",
        "Downregulated": "down"
    }

    kegg_df = get_kegg(
        kegg_map[kegg_regulation]
    )

    if kegg_df.empty:
        st.warning(
            "No KEGG enrichment results available for this selection."
        )

    else:
        st.markdown(
            f"### KEGG — {kegg_regulation}"
        )

        search_kegg = st.text_input(
            "Search pathway",
            placeholder="Example: focal adhesion",
            key="kegg_search"
        )

        kegg_table = kegg_df.copy()

        if search_kegg:
            kegg_table = kegg_table[
                kegg_table["Description"]
                .astype(str)
                .str.contains(
                    search_kegg,
                    case=False,
                    na=False
                )
            ]

        kegg_table = kegg_table.sort_values(
            "p.adjust",
            ascending=True
        )

        kegg_columns = [
            c for c in [
                "ID",
                "Description",
                "category",
                "subcategory",
                "GeneRatio",
                "Count",
                "pvalue",
                "p.adjust",
                "qvalue"
            ]
            if c in kegg_table.columns
        ]

        display_kegg = kegg_table[
            kegg_columns
        ].copy()

        st.dataframe(
            display_kegg.style.format({
                "pvalue": "{:.2e}",
                "p.adjust": "{:.2e}",
                "qvalue": "{:.2e}"
            }),
            width="stretch",
            hide_index=True,
            height=430
        )

        st.caption(
            f"{len(kegg_table)} KEGG pathways displayed."
        )

# ============================================================
# DOTPLOTS
# ============================================================
st.markdown(
    '<div class="section-title">📊 Enrichment Visualizations</div>',
    unsafe_allow_html=True
)

plot_left, plot_right = st.columns(
    2,
    gap="large"
)

with plot_left:

    st.markdown("### GO Biological Process")

    go_plot = (
        PLOT_DIR
        / "GO_BP_all_significant_dotplot.png"
    )

    if go_plot.exists():
        st.image(
            str(go_plot),
            width="stretch"
        )

    else:
        st.warning(
            "GO BP dotplot not found."
        )

with plot_right:

    st.markdown("### KEGG Pathways")

    kegg_plot = (
        PLOT_DIR
        / "KEGG_all_significant_dotplot.png"
    )

    if kegg_plot.exists():
        st.image(
            str(kegg_plot),
            width="stretch"
        )

    else:
        st.warning(
            "KEGG dotplot not found."
        )

# ============================================================
# TOP BIOLOGICAL SIGNALS
# ============================================================
st.markdown(
    '<div class="section-title">🏆 Strongest Enrichment Signals</div>',
    unsafe_allow_html=True
)

left, right = st.columns(
    2,
    gap="large"
)

with left:

    st.markdown("### Top GO Biological Processes")

    if not go_bp_all.empty:

        top_go = (
            go_bp_all
            .sort_values(
                "p.adjust",
                ascending=True
            )
            .head(10)
        )

        st.dataframe(
            top_go[
                [
                    "Description",
                    "Count",
                    "p.adjust"
                ]
            ].style.format({
                "p.adjust": "{:.2e}"
            }),
            width="stretch",
            hide_index=True
        )

with right:

    st.markdown("### Top KEGG Pathways")

    if not kegg_all.empty:

        top_kegg = (
            kegg_all
            .sort_values(
                "p.adjust",
                ascending=True
            )
            .head(10)
        )

        st.dataframe(
            top_kegg[
                [
                    "Description",
                    "Count",
                    "p.adjust"
                ]
            ].style.format({
                "p.adjust": "{:.2e}"
            }),
            width="stretch",
            hide_index=True
        )

# ============================================================
# BIOLOGICAL SUMMARY
# ============================================================
st.markdown(
    '<div class="section-title">🔬 Biological Interpretation</div>',
    unsafe_allow_html=True
)

st.markdown("""
<div class="info-box">

<strong>Functional interpretation of the RNAFlowX results</strong><br><br>

The enrichment results indicate coordinated transcriptional changes
across multiple biological systems rather than isolated changes in
individual genes.

GO Biological Process enrichment highlights processes associated with
vascular development, extracellular stimulus response and cellular
stress responses.

KEGG enrichment provides pathway-level context, including signaling,
cellular organization and stress-associated pathways.

These results complement the differential-expression analysis by
connecting individual DEGs to broader biological mechanisms.

</div>
""", unsafe_allow_html=True)

st.write("")

st.caption(
    "Functional enrichment: clusterProfiler • "
    "Significance criterion: adjusted p-value < 0.05"
)
