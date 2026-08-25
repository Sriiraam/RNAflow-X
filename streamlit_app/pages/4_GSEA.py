from utils.database import get_gsea

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go


# ============================================================
# PAGE CONFIG
# ============================================================

st.set_page_config(
    page_title="RNAFlowX | GSEA",
    layout="wide"
)


# ============================================================
# CSS
# ============================================================

st.markdown("""
<style>

.stApp{
    background:
        radial-gradient(circle at 8% 5%,
        rgba(124,58,237,.045), transparent 28%),
        #F8FAFC;
}

.gsea-hero{
    padding:42px 44px;
    border-radius:24px;

    background:
        linear-gradient(
            135deg,
            #031A46 0%,
            #173E8C 50%,
            #6D28D9 100%
        );

    color:white;

    box-shadow:
        0 18px 38px rgba(30,64,175,.18);

    margin-bottom:24px;
}

.gsea-hero h1{
    margin:0;
    font-size:46px;
    font-weight:850;
}

.gsea-hero h3{
    margin:8px 0 14px;
    color:#A7F3D0;
    font-size:20px;
}

.gsea-hero p{
    max-width:920px;
    color:#EEF2FF;
    font-size:15.5px;
    line-height:1.75;
}

.section-title{
    margin:36px 0 17px;
    font-size:29px;
    font-weight:850;
    color:#103C70;
}

.kpi{
    background:white;
    border:1px solid #E2E8F0;
    border-radius:18px;
    padding:20px;
    min-height:125px;
    box-shadow:0 6px 18px rgba(15,23,42,.055);
}

.kpi-label{
    font-size:13px;
    font-weight:700;
    color:#64748B;
}

.kpi-value{
    font-size:31px;
    font-weight:850;
    margin-top:8px;
}

.kpi-note{
    color:#94A3B8;
    font-size:12px;
    margin-top:7px;
}

.info-box{
    background:#F5F3FF;
    border-left:5px solid #7C3AED;
    padding:20px 22px;
    border-radius:14px;
    color:#334155;
    line-height:1.75;
}

</style>
""", unsafe_allow_html=True)


# ============================================================
# LOAD DATA
# ============================================================

try:
    gsea = get_gsea().copy()

except Exception as e:

    st.error(
        f"Unable to load GSEA results: {e}"
    )

    st.stop()


if gsea.empty:

    st.warning(
        "No GSEA results are available."
    )

    st.stop()


# ============================================================
# CLEAN DATA
# ============================================================

numeric_cols = [
    "setSize",
    "enrichmentScore",
    "NES",
    "pvalue",
    "p.adjust",
    "qvalue",
    "rank"
]

for col in numeric_cols:

    if col in gsea.columns:

        gsea[col] = pd.to_numeric(
            gsea[col],
            errors="coerce"
        )


gsea["padj_plot"] = (
    gsea["p.adjust"]
    .fillna(1)
    .clip(lower=1e-300)
)

gsea["neg_log10_padj"] = (
    -np.log10(
        gsea["padj_plot"]
    )
)


# ============================================================
# HERO
# ============================================================

st.markdown("""
<div class="gsea-hero">

<h1>Gene Set Enrichment Analysis</h1>

<h3>
Normalized Enrichment Score • GO Biological Process • Ranked Gene Sets
</h3>

<p>
Gene Set Enrichment Analysis identifies coordinated pathway-level
expression changes across a ranked gene list without requiring a strict
single-gene cutoff. Positive and negative normalized enrichment scores
highlight biological programs shifted toward opposite ends of the
differential-expression ranking.
</p>

</div>
""", unsafe_allow_html=True)


# ============================================================
# CONTROLS
# ============================================================

control_box = st.container(
    border=True
)

with control_box:

    c1, c2, c3 = st.columns(
        [1.3, 1.1, 1]
    )

    with c1:

        padj_threshold = st.slider(
            "Adjusted p-value threshold",
            min_value=0.001,
            max_value=0.10,
            value=0.05,
            step=0.001,
            format="%.3f"
        )

    with c2:

        top_n = st.slider(
            "Top pathways displayed",
            min_value=5,
            max_value=30,
            value=15,
            step=5
        )

    with c3:

        direction_filter = st.selectbox(
            "Direction",
            [
                "All Significant",
                "Positive NES",
                "Negative NES"
            ]
        )


# ============================================================
# FILTER
# ============================================================

significant = gsea[
    gsea["p.adjust"] < padj_threshold
].copy()

positive = significant[
    significant["NES"] > 0
].copy()

negative = significant[
    significant["NES"] < 0
].copy()


# ============================================================
# KPI
# ============================================================

st.markdown(
    '<div class="section-title">GSEA Overview</div>',
    unsafe_allow_html=True
)

k1, k2, k3, k4, k5 = st.columns(5)

cards = [
    (
        k1,
        "Gene Sets Tested",
        f"{len(gsea):,}",
        "#2563EB",
        "GO BP collections"
    ),
    (
        k2,
        "Significant Sets",
        f"{len(significant):,}",
        "#16A34A",
        f"padj < {padj_threshold:.3f}"
    ),
    (
        k3,
        "Positive NES",
        f"{len(positive):,}",
        "#DC2626",
        "Positive enrichment"
    ),
    (
        k4,
        "Negative NES",
        f"{len(negative):,}",
        "#2563EB",
        "Negative enrichment"
    ),
    (
        k5,
        "Strongest |NES|",
        (
            f"{significant['NES'].abs().max():.2f}"
            if not significant.empty
            else "0"
        ),
        "#7C3AED",
        "Largest pathway shift"
    )
]

for col, label, value, color, note in cards:

    with col:

        st.markdown(
            f"""
<div class="kpi" style="border-top:5px solid {color};">

<div class="kpi-label">
{label}
</div>

<div class="kpi-value"
style="color:{color};">
{value}
</div>

<div class="kpi-note">
{note}
</div>

</div>
""",
            unsafe_allow_html=True
        )


# ============================================================
# DIVERGING NES
# ============================================================

st.markdown(
    '<div class="section-title">Enrichment Direction Landscape</div>',
    unsafe_allow_html=True
)


top_positive = (
    positive
    .sort_values(
        "NES",
        ascending=False
    )
    .head(top_n)
)

top_negative = (
    negative
    .sort_values(
        "NES",
        ascending=True
    )
    .head(top_n)
)


diverging = pd.concat(
    [
        top_negative,
        top_positive
    ],
    ignore_index=True
)

if not diverging.empty:

    diverging["Direction"] = np.where(
        diverging["NES"] > 0,
        "Positive NES",
        "Negative NES"
    )

    diverging = diverging.sort_values(
        "NES"
    )


    nes_bar = px.bar(
        diverging,

        x="NES",
        y="Description",

        orientation="h",

        color="Direction",

        color_discrete_map={
            "Positive NES":
                "#EF4444",

            "Negative NES":
                "#2563EB"
        },

        custom_data=[
            "ID",
            "setSize",
            "p.adjust",
            "enrichmentScore",
            "rank"
        ],

        labels={
            "NES":
                "Normalized Enrichment Score",

            "Description":
                ""
        }
    )


    nes_bar.update_traces(
        hovertemplate=
        "<b>%{y}</b>"
        "<br>ID: %{customdata[0]}"
        "<br>NES: %{x:.3f}"
        "<br>Set size: %{customdata[1]}"
        "<br>padj: %{customdata[2]:.3e}"
        "<br>Enrichment score: %{customdata[3]:.3f}"
        "<br>Rank: %{customdata[4]}"
        "<extra></extra>"
    )


    nes_bar.add_vline(
        x=0,
        line_dash="dot",
        line_color="#64748B"
    )


    nes_bar.update_layout(
        template="plotly_white",
        height=720,
        legend_title_text="",
        margin=dict(
            l=10,
            r=20,
            t=20,
            b=20
        )
    )


    st.plotly_chart(
        nes_bar,
        width="stretch"
    )


# ============================================================
# DISTRIBUTIONS
# ============================================================

st.markdown(
    '<div class="section-title">NES Distribution</div>',
    unsafe_allow_html=True
)

d1, d2 = st.columns(
    2,
    gap="large"
)


# ============================================================
# HISTOGRAM + DENSITY
# ============================================================

with d1:

    st.markdown(
        "### NES Histogram"
    )


    hist = px.histogram(
        significant,

        x="NES",

        nbins=30,

        color=np.where(
            significant["NES"] >= 0,
            "Positive",
            "Negative"
        ),

        color_discrete_map={
            "Positive":
                "#EF4444",

            "Negative":
                "#2563EB"
        },

        marginal="rug",

        labels={
            "NES":
                "Normalized Enrichment Score"
        }
    )


    hist.update_layout(
        template="plotly_white",
        height=470,
        showlegend=True,
        legend_title_text=""
    )


    st.plotly_chart(
        hist,
        width="stretch"
    )


# ============================================================
# KDE-LIKE DENSITY
# ============================================================

with d2:

    st.markdown(
        "### Density View"
    )


    density = go.Figure()


    for label, subset, color in [

        (
            "Positive NES",
            positive,
            "#EF4444"
        ),

        (
            "Negative NES",
            negative,
            "#2563EB"
        )

    ]:

        if not subset.empty:

            density.add_trace(
                go.Histogram(
                    x=subset["NES"],

                    name=label,

                    histnorm="probability density",

                    nbinsx=25,

                    opacity=.45,

                    marker_color=color
                )
            )


    density.update_layout(
        barmode="overlay",

        template="plotly_white",

        height=470,

        xaxis_title=
            "Normalized Enrichment Score",

        yaxis_title=
            "Density",

        legend_title_text=""
    )


    st.plotly_chart(
        density,
        width="stretch"
    )


# ============================================================
# NES VS SIGNIFICANCE
# ============================================================

st.markdown(
    '<div class="section-title">NES vs Statistical Significance</div>',
    unsafe_allow_html=True
)


scatter = px.scatter(
    significant,

    x="NES",

    y="neg_log10_padj",

    size="setSize",

    color=np.where(
        significant["NES"] > 0,
        "Positive NES",
        "Negative NES"
    ),

    color_discrete_map={
        "Positive NES":
            "#EF4444",

        "Negative NES":
            "#2563EB"
    },

    hover_name="Description",

    hover_data={
        "ID": True,
        "setSize": True,
        "p.adjust": ":.2e",
        "NES": ":.3f",
        "enrichmentScore": ":.3f",
        "neg_log10_padj": False
    },

    labels={
        "NES":
            "Normalized Enrichment Score",

        "neg_log10_padj":
            "-log10 adjusted p-value"
    }
)


scatter.add_vline(
    x=0,
    line_dash="dot",
    line_color="#64748B"
)


scatter.update_layout(
    template="plotly_white",
    height=550,
    legend_title_text=""
)


st.plotly_chart(
    scatter,
    width="stretch"
)


# ============================================================
# RUG / STRIP STYLE
# ============================================================

st.markdown(
    '<div class="section-title">Pathway NES Strip</div>',
    unsafe_allow_html=True
)


strip_df = significant.copy()

strip_df["Direction"] = np.where(
    strip_df["NES"] > 0,
    "Positive NES",
    "Negative NES"
)


strip = px.strip(
    strip_df,

    x="NES",

    y="Direction",

    color="Direction",

    hover_name="Description",

    hover_data={
        "setSize": True,
        "p.adjust": ":.2e",
        "NES": ":.3f"
    },

    color_discrete_map={
        "Positive NES":
            "#EF4444",

        "Negative NES":
            "#2563EB"
    }
)


strip.update_traces(
    marker=dict(
        size=8,
        opacity=.55
    )
)


strip.add_vline(
    x=0,
    line_dash="dot",
    line_color="#64748B"
)


strip.update_layout(
    template="plotly_white",
    height=330,
    showlegend=False
)


st.plotly_chart(
    strip,
    width="stretch"
)


# ============================================================
# GENE SET EXPLORER
# ============================================================

st.markdown(
    '<div class="section-title">GSEA Gene-Set Explorer</div>',
    unsafe_allow_html=True
)


f1, f2 = st.columns(
    [2, 1]
)


with f1:

    search = st.text_input(
        "Search biological process",
        placeholder="Example: cholesterol biosynthetic process"
    )


with f2:

    explorer_direction = st.selectbox(
        "Explorer direction",
        [
            "All Significant",
            "Positive NES",
            "Negative NES"
        ],
        key="gsea_explorer_direction"
    )


if explorer_direction == "Positive NES":

    table = positive.copy()

elif explorer_direction == "Negative NES":

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
    "p.adjust"
)


display_columns = [
    "ID",
    "Description",
    "setSize",
    "enrichmentScore",
    "NES",
    "pvalue",
    "p.adjust",
    "qvalue",
    "rank"
]


st.dataframe(
    table[
        display_columns
    ],

    width="stretch",
    hide_index=True,
    height=430,

    column_config={

        "NES":
            st.column_config.NumberColumn(
                format="%.2f"
            ),

        "enrichmentScore":
            st.column_config.NumberColumn(
                format="%.3f"
            ),

        "pvalue":
            st.column_config.NumberColumn(
                format="%.2e"
            ),

        "p.adjust":
            st.column_config.NumberColumn(
                format="%.2e"
            ),

        "qvalue":
            st.column_config.NumberColumn(
                format="%.2e"
            )
    }
)


st.caption(
    f"{len(table):,} gene sets match "
    "the current filters."
)


# ============================================================
# TOP POSITIVE / NEGATIVE
# ============================================================

st.markdown(
    '<div class="section-title">Strongest Enrichment Signals</div>',
    unsafe_allow_html=True
)


t1, t2 = st.columns(
    2,
    gap="large"
)


with t1:

    st.markdown(
        "### 🔴 Top Positive NES"
    )

    st.dataframe(
        top_positive[
            [
                "Description",
                "setSize",
                "NES",
                "p.adjust"
            ]
        ],

        width="stretch",
        hide_index=True,

        column_config={
            "NES":
                st.column_config.NumberColumn(
                    format="%.2f"
                ),

            "p.adjust":
                st.column_config.NumberColumn(
                    format="%.2e"
                )
        }
    )


with t2:

    st.markdown(
        "### 🔵 Top Negative NES"
    )

    st.dataframe(
        top_negative[
            [
                "Description",
                "setSize",
                "NES",
                "p.adjust"
            ]
        ],

        width="stretch",
        hide_index=True,

        column_config={
            "NES":
                st.column_config.NumberColumn(
                    format="%.2f"
                ),

            "p.adjust":
                st.column_config.NumberColumn(
                    format="%.2e"
                )
        }
    )


# ============================================================
# LEADING EDGE EXPLORER
# ============================================================

st.markdown(
    '<div class="section-title">Leading-Edge Explorer</div>',
    unsafe_allow_html=True
)


available_sets = (
    significant
    .sort_values(
        "NES",
        key=lambda s: s.abs(),
        ascending=False
    )
    ["Description"]
    .tolist()
)


if available_sets:

    selected_pathway = st.selectbox(
        "Select pathway",
        available_sets
    )


    selected = significant[
        significant["Description"]
        == selected_pathway
    ].iloc[0]


    c1, c2, c3, c4 = st.columns(4)

    c1.metric(
        "NES",
        f"{selected['NES']:.2f}"
    )

    c2.metric(
        "Set Size",
        f"{int(selected['setSize'])}"
    )

    c3.metric(
        "Adjusted p-value",
        f"{selected['p.adjust']:.2e}"
    )

    c4.metric(
        "Rank",
        f"{int(selected['rank'])}"
    )


    st.markdown(
        "#### Leading-edge summary"
    )

    st.info(
        str(
            selected[
                "leading_edge"
            ]
        )
    )


    st.markdown(
        "#### Core enrichment genes"
    )


    core_genes = str(
        selected[
            "core_enrichment"
        ]
    ).split("/")


    core_gene_df = pd.DataFrame(
        {
            "Core Enrichment Gene":
                core_genes
        }
    )


    st.dataframe(
        core_gene_df,
        width="stretch",
        hide_index=True,
        height=300
    )


# ============================================================
# INTERPRETATION
# ============================================================

st.markdown(
    '<div class="section-title">Biological Interpretation</div>',
    unsafe_allow_html=True
)


strongest_positive = (
    positive
    .sort_values(
        "NES",
        ascending=False
    )
    .iloc[0]
    if not positive.empty
    else None
)


strongest_negative = (
    negative
    .sort_values(
        "NES",
        ascending=True
    )
    .iloc[0]
    if not negative.empty
    else None
)


positive_text = (
    strongest_positive["Description"]
    if strongest_positive is not None
    else "No significant positive pathway"
)

negative_text = (
    strongest_negative["Description"]
    if strongest_negative is not None
    else "No significant negative pathway"
)


st.markdown(
    f"""
<div class="info-box">

<strong>Strongest positive enrichment:</strong>
{positive_text}

<br><br>

<strong>Strongest negative enrichment:</strong>
{negative_text}

<br><br>

Positive NES values indicate pathways concentrated toward
the upregulated end of the ranked gene list. Negative NES
values indicate pathways concentrated toward the
downregulated end.

GSEA complements conventional over-representation analysis
because it detects coordinated pathway-level shifts even
when individual genes do not pass a strict differential
expression threshold.

</div>
""",
    unsafe_allow_html=True
)


# ============================================================
# DOWNLOAD
# ============================================================

st.markdown(
    '<div class="section-title">Download GSEA Results</div>',
    unsafe_allow_html=True
)


download_csv = (
    significant
    .to_csv(
        index=False
    )
    .encode("utf-8")
)


st.download_button(
    "Download Significant GSEA Results",
    data=download_csv,
    file_name="gsea_significant_current_threshold.csv",
    mime="text/csv",
    width="stretch"
)


st.write("")

st.caption(
    "GSEA source: clusterProfiler • "
    "GO Biological Process • "
    "Multiple-testing correction: Benjamini-Hochberg FDR"
)
