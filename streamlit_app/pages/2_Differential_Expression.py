from utils.database import get_differential_expression

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from pathlib import Path
from sklearn.decomposition import PCA


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

APP_ROOT = Path(__file__).resolve().parents[1]

DOWNLOAD_DIR = APP_ROOT / "data" / "downloads"

NORMALIZED_COUNTS_FILE = (
    DOWNLOAD_DIR / "normalized_counts.csv"
)

DE_FILE = (
    DOWNLOAD_DIR / "differential_expression.csv"
)


# ============================================================
# CSS
# ============================================================

st.markdown(
    """
<style>

.stApp{
    background:
        radial-gradient(circle at 8% 5%,
        rgba(37,99,235,.045), transparent 28%),
        #F8FAFC;
}

/* HERO */

.de-hero{
    padding:40px 44px;
    border-radius:24px;

    background:
        linear-gradient(
            135deg,
            #031A46 0%,
            #063D91 52%,
            #075FCB 100%
        );

    color:white;

    box-shadow:
        0 18px 38px rgba(15,74,180,.18);

    margin-bottom:22px;
}

.de-hero h1{
    margin:0;
    font-size:46px;
    font-weight:850;
    letter-spacing:-.6px;
}

.de-hero h3{
    color:#55F1DB;
    font-size:20px;
    margin:8px 0 14px;
}

.de-hero p{
    color:#E7F1FF;
    max-width:900px;
    line-height:1.75;
    font-size:15.5px;
}


/* SECTION */

.de-section{
    margin:36px 0 17px;
    font-size:29px;
    font-weight:850;
    color:#103C70;
}


/* KPI */

.de-kpi{
    background:white;

    border:1px solid #E2E8F0;
    border-radius:18px;

    padding:20px;

    min-height:132px;

    box-shadow:
        0 6px 18px rgba(15,23,42,.055);
}

.de-kpi-label{
    font-size:13px;
    font-weight:700;
    color:#64748B;
}

.de-kpi-value{
    font-size:31px;
    line-height:1.1;
    font-weight:850;
    margin-top:8px;
}

.de-kpi-note{
    font-size:12px;
    margin-top:7px;
    color:#94A3B8;
}


/* INFO */

.de-info{
    background:white;

    border:1px solid #E2E8F0;
    border-radius:17px;

    padding:20px;

    line-height:1.7;

    box-shadow:
        0 5px 16px rgba(15,23,42,.045);
}

.de-info h3{
    color:#123D72;
    margin-top:0;
}


/* MINI LABEL */

.de-caption{
    margin-top:-7px;
    text-align:center;
    font-size:12px;
    color:#64748B;
}


/* STAT BOX */

.stat-box{
    background:#EFF6FF;
    border-left:5px solid #2563EB;
    border-radius:14px;
    padding:20px 22px;
    color:#334155;
    line-height:1.75;
}

</style>
""",
    unsafe_allow_html=True
)


# ============================================================
# HERO
# ============================================================

st.markdown(
    """
<div class="de-hero">

<h1>Differential Expression Analysis</h1>

<h3>
DESeq2 • Statistical Transcriptomics • Interactive Exploration
</h3>

<p>
Explore transcriptional differences between control and PFOS-treated
samples using DESeq2 results. Interactive significance thresholds
automatically update the differential-expression classification,
volcano plot, MA plot, summary statistics and gene explorer.
</p>

</div>
""",
    unsafe_allow_html=True
)


# ============================================================
# LOAD DATA
# ============================================================

if not DE_FILE.exists():

    st.error(
        "Packaged differential-expression results were not found."
    )

    st.stop()


if not NORMALIZED_COUNTS_FILE.exists():

    st.error(
        "Packaged normalized-count matrix was not found."
    )

    st.stop()


deg = get_differential_expression().copy()

counts = pd.read_csv(
    NORMALIZED_COUNTS_FILE
)


# ============================================================
# VALIDATE DE DATA
# ============================================================

if "gene_id" in deg.columns and "Gene" not in deg.columns:
    deg = deg.rename(columns={"gene_id": "Gene"})

required = {
    "Gene",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
}

missing = required - set(deg.columns)

if missing:

    st.error(
        "Missing DESeq2 columns: "
        + ", ".join(sorted(missing))
    )

    st.stop()


for col in [
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
]:

    deg[col] = pd.to_numeric(
        deg[col],
        errors="coerce"
    )


deg["padj_plot"] = (
    deg["padj"]
    .fillna(1)
    .clip(lower=1e-300)
)

deg["neg_log10_padj"] = (
    -np.log10(deg["padj_plot"])
)


# ============================================================
# NORMALIZED COUNTS
# ============================================================

gene_column = counts.columns[0]

counts = counts.rename(
    columns={
        gene_column: "Gene"
    }
)

counts["Gene"] = counts["Gene"].astype(str)

counts = counts.set_index("Gene")

sample_columns = counts.columns.tolist()

counts = counts.apply(
    pd.to_numeric,
    errors="coerce"
).fillna(0)


# ============================================================
# SAMPLE METADATA
# ============================================================

sample_metadata = []

for sample in sample_columns:

    lower = sample.lower()

    if "control" in lower:

        condition = "Control"

    elif "pfos" in lower:

        condition = "PFOS 50 µM"

    else:

        condition = "Other"

    sample_metadata.append(
        {
            "Sample": sample,
            "Condition": condition
        }
    )

sample_metadata = pd.DataFrame(
    sample_metadata
)


# ============================================================
# THRESHOLDS
# ============================================================

control_box = st.container(border=True)

with control_box:

    c1, c2, c3 = st.columns(
        [1.4, 1.4, 0.8],
        gap="large"
    )

    with c1:

        st.markdown("#### 🎛️ Statistical Significance")

        padj_threshold = st.slider(
            "Adjusted p-value (padj)",
            min_value=0.001,
            max_value=0.10,
            value=0.05,
            step=0.001,
            format="%.3f"
        )

        st.caption(
            "Benjamini–Hochberg FDR threshold"
        )

    with c2:

        st.markdown("#### 📐 Effect Size")

        lfc_threshold = st.slider(
            "|log2 Fold Change|",
            min_value=0.0,
            max_value=5.0,
            value=1.0,
            step=0.1
        )

        st.caption(
            "Minimum absolute transcriptional change"
        )

    with c3:

        st.markdown("#### 🧪 Comparison")

        st.markdown(
            """
            **PFOS 50 µM**  
            vs  
            **Control**
            """
        )

        st.caption(
            "DESeq2 • BH-FDR"
        )

st.caption(
    "🟢 Changing either threshold automatically updates "
    "the DEG counts, Volcano plot, MA plot and gene tables."
)

# ============================================================
# CLASSIFY GENES
# ============================================================

deg["Regulation"] = "Not significant"

up_mask = (
    (deg["padj"] < padj_threshold)
    &
    (deg["log2FoldChange"] >= lfc_threshold)
)

down_mask = (
    (deg["padj"] < padj_threshold)
    &
    (deg["log2FoldChange"] <= -lfc_threshold)
)

deg.loc[
    up_mask,
    "Regulation"
] = "Upregulated"

deg.loc[
    down_mask,
    "Regulation"
] = "Downregulated"


significant = deg[
    deg["Regulation"] != "Not significant"
].copy()

upregulated = deg[
    deg["Regulation"] == "Upregulated"
].copy()

downregulated = deg[
    deg["Regulation"] == "Downregulated"
].copy()


# ============================================================
# KPI CARDS
# ============================================================

st.markdown(
    '<div class="de-section">Differential Expression Overview</div>',
    unsafe_allow_html=True
)


k1, k2, k3, k4, k5 = st.columns(5)


cards = [

    (
        k1,
        "Genes Tested",
        f"{len(deg):,}",
        "#2563EB",
        "DESeq2 results"
    ),

    (
        k2,
        "Significant DEGs",
        f"{len(significant):,}",
        "#16A34A",
        f"padj < {padj_threshold:.3f}"
    ),

    (
        k3,
        "Upregulated",
        f"{len(upregulated):,}",
        "#DC2626",
        "Higher in PFOS"
    ),

    (
        k4,
        "Downregulated",
        f"{len(downregulated):,}",
        "#2563EB",
        "Lower in PFOS"
    ),

    (
        k5,
        "Samples",
        f"{len(sample_columns)}",
        "#7C3AED",
        "2 control • 2 PFOS"
    )

]


for column, label, value, color, note in cards:

    with column:

        st.markdown(
            f"""
<div class="de-kpi"
style="border-top:5px solid {color};">

<div class="de-kpi-label">
{label}
</div>

<div class="de-kpi-value"
style="color:{color};">
{value}
</div>

<div class="de-kpi-note">
{note}
</div>

</div>
""",
            unsafe_allow_html=True
        )


# ============================================================
# PCA
# ============================================================

log_counts = np.log2(
    counts + 1
).T


pca_model = PCA(
    n_components=2
)

pca_coordinates = (
    pca_model.fit_transform(
        log_counts
    )
)


pca_df = pd.DataFrame(
    {
        "Sample":
            log_counts.index,

        "PC1":
            pca_coordinates[:, 0],

        "PC2":
            pca_coordinates[:, 1]
    }
)


pca_df = pca_df.merge(
    sample_metadata,
    on="Sample",
    how="left"
)


pc1_variance = (
    pca_model
    .explained_variance_ratio_[0]
    * 100
)

pc2_variance = (
    pca_model
    .explained_variance_ratio_[1]
    * 100
)


# ============================================================
# EXPRESSION LANDSCAPE
# ============================================================

st.markdown(
    '<div class="de-section">Expression Landscape</div>',
    unsafe_allow_html=True
)


pca_col, volcano_col, ma_col = st.columns(
    3,
    gap="large"
)


# ============================================================
# PCA PLOT
# ============================================================

with pca_col:

    st.markdown("### Principal Component Analysis")

    pca_fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color="Condition",
        text="Sample",

        color_discrete_map={
            "Control": "#2563EB",
            "PFOS 50 µM": "#EF4444"
        },

        hover_data={
            "Sample": True,
            "Condition": True,
            "PC1": ":.2f",
            "PC2": ":.2f"
        },

        labels={
            "PC1":
                f"PC1 ({pc1_variance:.1f}%)",

            "PC2":
                f"PC2 ({pc2_variance:.1f}%)"
        }
    )


    pca_fig.update_traces(
        marker=dict(
            size=15,
            line=dict(
                width=1.5,
                color="white"
            )
        ),
        textposition="top center"
    )


    pca_fig.update_layout(
        template="plotly_white",
        height=500,
        legend_title_text="",
        margin=dict(
            l=25,
            r=15,
            t=25,
            b=20
        )
    )


    pca_fig.add_hline(
        y=0,
        line_dash="dot",
        line_color="#94A3B8"
    )

    pca_fig.add_vline(
        x=0,
        line_dash="dot",
        line_color="#94A3B8"
    )


    st.plotly_chart(
        pca_fig,
        width="stretch"
    )

    st.markdown(
        '<div class="de-caption">'
        'PCA derived from log2-normalized expression counts'
        '</div>',
        unsafe_allow_html=True
    )


# ============================================================
# VOLCANO
# ============================================================

with volcano_col:

    st.markdown("### Volcano Plot")

    volcano_fig = px.scatter(
        deg,

        x="log2FoldChange",
        y="neg_log10_padj",

        color="Regulation",

        color_discrete_map={
            "Upregulated":
                "#FF304F",

            "Downregulated":
                "#2563FF",

            "Not significant":
                "#CBD5E1"
        },

        category_orders={
            "Regulation": [
                "Not significant",
                "Downregulated",
                "Upregulated"
            ]
        },

        custom_data=[
            "Gene",
            "baseMean",
            "padj",
            "pvalue",
            "Regulation"
        ],

        labels={
            "log2FoldChange":
                "log2 Fold Change",

            "neg_log10_padj":
                "-log10 adjusted p-value"
        }
    )


    volcano_fig.update_traces(
        marker=dict(
            size=6,
            opacity=.72
        ),

        hovertemplate=
        "<b>%{customdata[0]}</b>"
        "<br>Regulation: %{customdata[4]}"
        "<br>log2FC: %{x:.3f}"
        "<br>-log10(padj): %{y:.3f}"
        "<br>baseMean: %{customdata[1]:.2f}"
        "<br>padj: %{customdata[2]:.3e}"
        "<br>p-value: %{customdata[3]:.3e}"
        "<extra></extra>"
    )


    volcano_fig.add_vline(
        x=lfc_threshold,
        line_dash="dash",
        line_color="#EF4444",
        opacity=.65
    )

    volcano_fig.add_vline(
        x=-lfc_threshold,
        line_dash="dash",
        line_color="#2563EB",
        opacity=.65
    )


    volcano_fig.add_hline(
        y=-np.log10(
            padj_threshold
        ),
        line_dash="dash",
        line_color="#64748B",
        opacity=.7
    )


    volcano_fig.update_layout(
        template="plotly_white",
        height=500,
        legend_title_text="",
        margin=dict(
            l=25,
            r=15,
            t=25,
            b=20
        )
    )


    st.plotly_chart(
        volcano_fig,
        width="stretch"
    )


    st.markdown(
        '<div class="de-caption">'
        'Red = upregulated • Blue = downregulated • Grey = not significant'
        '</div>',
        unsafe_allow_html=True
    )


# ============================================================
# MA PLOT
# ============================================================

with ma_col:

    st.markdown("### MA Plot")

    ma_data = deg[
        (deg["baseMean"] > 0)
        &
        deg["baseMean"].notna()
        &
        deg["log2FoldChange"].notna()
    ].copy()


    ma_fig = px.scatter(
        ma_data,

        x="baseMean",
        y="log2FoldChange",

        color="Regulation",

        color_discrete_map={
            "Upregulated":
                "#FF304F",

            "Downregulated":
                "#2563FF",

            "Not significant":
                "#CBD5E1"
        },

        custom_data=[
            "Gene",
            "padj",
            "baseMean",
            "Regulation"
        ],

        log_x=True,

        labels={
            "baseMean":
                "Mean normalized counts",

            "log2FoldChange":
                "log2 Fold Change"
        }
    )


    ma_fig.update_traces(
        marker=dict(
            size=5,
            opacity=.65
        ),

        hovertemplate=
        "<b>%{customdata[0]}</b>"
        "<br>Regulation: %{customdata[3]}"
        "<br>baseMean: %{customdata[2]:.2f}"
        "<br>log2FC: %{y:.3f}"
        "<br>padj: %{customdata[1]:.3e}"
        "<extra></extra>"
    )


    ma_fig.add_hline(
        y=0,
        line_dash="dot",
        line_color="#64748B"
    )

    ma_fig.add_hline(
        y=lfc_threshold,
        line_dash="dash",
        line_color="#EF4444",
        opacity=.55
    )

    ma_fig.add_hline(
        y=-lfc_threshold,
        line_dash="dash",
        line_color="#2563EB",
        opacity=.55
    )


    ma_fig.update_layout(
        template="plotly_white",
        height=500,
        legend_title_text="",
        margin=dict(
            l=25,
            r=15,
            t=25,
            b=20
        )
    )


    st.plotly_chart(
        ma_fig,
        width="stretch"
    )


    st.markdown(
        '<div class="de-caption">'
        'Mean expression versus estimated transcriptional effect'
        '</div>',
        unsafe_allow_html=True
    )


# ============================================================
# EFFECT DISTRIBUTION
# ============================================================

st.markdown(
    '<div class="de-section">Expression Effect Distribution</div>',
    unsafe_allow_html=True
)


dist_left, dist_right = st.columns(
    2,
    gap="large"
)


# ============================================================
# SIGNIFICANT EFFECT DISTRIBUTION
# ============================================================

with dist_left:

    st.markdown(
        "### Significant Differential Signals"
    )


    if not significant.empty:

        distribution_fig = go.Figure()


        for regulation, color in [

            (
                "Upregulated",
                "#FF4D5A"
            ),

            (
                "Downregulated",
                "#2878FF"
            )

        ]:

            subset = significant[
                significant["Regulation"]
                == regulation
            ]


            distribution_fig.add_trace(
                go.Violin(
                    y=subset[
                        "log2FoldChange"
                    ],

                    name=regulation,

                    box_visible=True,

                    meanline_visible=True,

                    points="all",

                    jitter=.32,

                    pointpos=0,

                    marker=dict(
                        size=4,
                        opacity=.35,
                        color=color
                    ),

                    line_color=color,

                    fillcolor=color,

                    opacity=.55
                )
            )


        distribution_fig.update_layout(
            template="plotly_white",
            height=450,
            yaxis_title="log2 Fold Change",
            xaxis_title="",
            showlegend=False,
            margin=dict(
                l=25,
                r=15,
                t=25,
                b=20
            )
        )


        st.plotly_chart(
            distribution_fig,
            width="stretch"
        )


    else:

        st.info(
            "No genes satisfy the current "
            "significance thresholds."
        )


# ============================================================
# ALL GENE EFFECT DISTRIBUTION
# ============================================================

with dist_right:

    st.markdown(
        "### Fold-Change Distribution"
    )


    valid_effects = deg[
        "log2FoldChange"
    ].dropna()


    effect_fig = go.Figure()


    effect_fig.add_trace(
        go.Violin(
            y=valid_effects,

            name="All genes",

            box_visible=True,

            meanline_visible=True,

            points="all",

            jitter=.35,

            pointpos=0,

            marker=dict(
                size=3,
                opacity=.18,
                color="#8B5CF6"
            ),

            line_color="#7C3AED",

            fillcolor="#A78BFA",

            opacity=.55
        )
    )


    effect_fig.add_hline(
        y=0,
        line_dash="dot",
        line_color="#64748B"
    )


    effect_fig.update_layout(
        template="plotly_white",
        height=450,
        yaxis_title="log2 Fold Change",
        showlegend=False,
        margin=dict(
            l=25,
            r=15,
            t=25,
            b=20
        )
    )


    st.plotly_chart(
        effect_fig,
        width="stretch"
    )


# ============================================================
# GENE EXPLORER
# ============================================================

st.markdown(
    '<div class="de-section">Significant Gene Explorer</div>',
    unsafe_allow_html=True
)


f1, f2, f3 = st.columns(
    [1.5, 1, 1]
)


with f1:

    gene_search = st.text_input(
        "Search Gene ID",
        placeholder="Example: ENSG000001..."
    )


with f2:

    regulation_filter = st.selectbox(
        "Regulation",
        [
            "All significant",
            "Upregulated",
            "Downregulated"
        ]
    )


with f3:

    sort_option = st.selectbox(
        "Sort by",
        [
            "padj",
            "Absolute log2FC",
            "baseMean"
        ]
    )


if regulation_filter == "Upregulated":

    gene_table = (
        upregulated.copy()
    )

elif regulation_filter == "Downregulated":

    gene_table = (
        downregulated.copy()
    )

else:

    gene_table = (
        significant.copy()
    )


if gene_search:

    gene_table = gene_table[
        gene_table["Gene"]
        .astype(str)
        .str.contains(
            gene_search,
            case=False,
            na=False
        )
    ]


if sort_option == "Absolute log2FC":

    gene_table["abs_lfc"] = (
        gene_table[
            "log2FoldChange"
        ].abs()
    )

    gene_table = gene_table.sort_values(
        "abs_lfc",
        ascending=False
    )

elif sort_option == "baseMean":

    gene_table = gene_table.sort_values(
        "baseMean",
        ascending=False
    )

else:

    gene_table = gene_table.sort_values(
        "padj",
        ascending=True,
        na_position="last"
    )


display_columns = [
    "Gene",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj",
    "Regulation"
]


st.dataframe(
    gene_table[
        display_columns
    ],
    width="stretch",
    hide_index=True,
    height=430,

    column_config={

        "baseMean":
            st.column_config.NumberColumn(
                format="%.1f"
            ),

        "log2FoldChange":
            st.column_config.NumberColumn(
                format="%.2f"
            ),

        "pvalue":
            st.column_config.NumberColumn(
                format="%.2e"
            ),

        "padj":
            st.column_config.NumberColumn(
                format="%.2e"
            )
    }
)


st.caption(
    f"{len(gene_table):,} genes match "
    "the current filters."
)


# ============================================================
# TOP GENES
# ============================================================

st.markdown(
    '<div class="de-section">Strongest Differential Signals</div>',
    unsafe_allow_html=True
)


top_up_col, top_down_col = st.columns(
    2,
    gap="large"
)


with top_up_col:

    st.markdown(
        "### 🔴 Top Upregulated Genes"
    )

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
        ],

        width="stretch",
        hide_index=True,

        column_config={

            "log2FoldChange":
                st.column_config.NumberColumn(
                    format="%.2f"
                ),

            "padj":
                st.column_config.NumberColumn(
                    format="%.2e"
                )
        }
    )


with top_down_col:

    st.markdown(
        "### 🔵 Top Downregulated Genes"
    )

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
        ],

        width="stretch",
        hide_index=True,

        column_config={

            "log2FoldChange":
                st.column_config.NumberColumn(
                    format="%.2f"
                ),

            "padj":
                st.column_config.NumberColumn(
                    format="%.2e"
                )
        }
    )


# ============================================================
# GENE EXPRESSION VIEWER
# ============================================================

st.markdown(
    '<div class="de-section">Gene Expression Viewer</div>',
    unsafe_allow_html=True
)


available_genes = (
    significant["Gene"]
    .astype(str)
    .tolist()
)


# Match only genes present in normalized counts
available_genes = [
    gene
    for gene in available_genes
    if gene in counts.index
]


if not available_genes:

    st.info(
        "No significant genes could be matched "
        "to the normalized-count matrix."
    )

else:

    selected_gene = st.selectbox(
        "Select significant gene",
        available_genes
    )


    selected_values = (
        counts.loc[
            selected_gene
        ]
        .reset_index()
    )


    selected_values.columns = [
        "Sample",
        "Normalized Count"
    ]


    selected_values = (
        selected_values.merge(
            sample_metadata,
            on="Sample",
            how="left"
        )
    )


    selected_values[
        "log2 Normalized Count"
    ] = np.log2(
        selected_values[
            "Normalized Count"
        ] + 1
    )


    selected_stats = deg[
        deg["Gene"]
        .astype(str)
        == selected_gene
    ]


    if not selected_stats.empty:

        gene_row = (
            selected_stats.iloc[0]
        )

        m1, m2, m3, m4 = st.columns(4)

        m1.metric(
            "Base Mean",
            f"{gene_row['baseMean']:,.1f}"
        )

        m2.metric(
            "log2FC",
            f"{gene_row['log2FoldChange']:.2f}"
        )

        m3.metric(
            "Adjusted p-value",
            f"{gene_row['padj']:.2e}"
        )

        m4.metric(
            "Regulation",
            gene_row["Regulation"]
        )


    viewer1, viewer2, viewer3 = (
        st.columns(
            3,
            gap="large"
        )
    )


    # VIOLIN
    with viewer1:

        st.markdown("#### Violin + Points")

        gene_violin = px.violin(
            selected_values,

            x="Condition",

            y="log2 Normalized Count",

            color="Condition",

            box=True,

            points="all",

            hover_data=[
                "Sample",
                "Normalized Count"
            ],

            color_discrete_map={
                "Control":
                    "#2563EB",

                "PFOS 50 µM":
                    "#EF4444"
            }
        )


        gene_violin.update_layout(
            template="plotly_white",
            height=420,
            showlegend=False,
            margin=dict(
                l=15,
                r=10,
                t=20,
                b=20
            )
        )


        st.plotly_chart(
            gene_violin,
            width="stretch"
        )


    # BOX + STRIP
    with viewer2:

        st.markdown("#### Box + Individual Samples")

        gene_box = px.box(
            selected_values,

            x="Condition",

            y="log2 Normalized Count",

            color="Condition",

            points="all",

            hover_data=[
                "Sample",
                "Normalized Count"
            ],

            color_discrete_map={
                "Control":
                    "#2563EB",

                "PFOS 50 µM":
                    "#EF4444"
            }
        )


        gene_box.update_traces(
            jitter=.32
        )


        gene_box.update_layout(
            template="plotly_white",
            height=420,
            showlegend=False,
            margin=dict(
                l=15,
                r=10,
                t=20,
                b=20
            )
        )


        st.plotly_chart(
            gene_box,
            width="stretch"
        )


    # STRIP
    with viewer3:

        st.markdown("#### Sample-Level Strip Plot")

        strip_fig = px.strip(
            selected_values,

            x="Condition",

            y="log2 Normalized Count",

            color="Condition",

            hover_name="Sample",

            hover_data={
                "Normalized Count":
                    ":.2f"
            },

            color_discrete_map={
                "Control":
                    "#2563EB",

                "PFOS 50 µM":
                    "#EF4444"
            }
        )


        strip_fig.update_traces(
            marker=dict(
                size=13,
                opacity=.85,
                line=dict(
                    width=1,
                    color="white"
                )
            )
        )


        strip_fig.update_layout(
            template="plotly_white",
            height=420,
            showlegend=False,
            margin=dict(
                l=15,
                r=10,
                t=20,
                b=20
            )
        )


        st.plotly_chart(
            strip_fig,
            width="stretch"
        )


    st.caption(
        "Only two biological replicates are available "
        "per condition; gene-level distribution plots "
        "should therefore be interpreted descriptively."
    )


# ============================================================
# STATISTICAL INTERPRETATION
# ============================================================

st.markdown(
    '<div class="de-section">Statistical Interpretation</div>',
    unsafe_allow_html=True
)


st.markdown(
    f"""
<div class="stat-box">

The DESeq2 analysis evaluated
<strong>{len(deg):,}</strong> genes.

Using:

<strong>padj &lt; {padj_threshold:.3f}</strong>

and

<strong>|log2FC| ≥ {lfc_threshold:.1f}</strong>

<br><br>

🔴 <strong>{len(upregulated):,}</strong>
genes are upregulated.

<br>

🔵 <strong>{len(downregulated):,}</strong>
genes are downregulated.

<br>

🧬 <strong>{len(significant):,}</strong>
genes satisfy the complete differential-expression definition.

<br><br>

These significant genes provide the primary input for downstream
GO, KEGG and gene-set enrichment analyses.

</div>
""",
    unsafe_allow_html=True
)


# ============================================================
# DOWNLOADS
# ============================================================

st.markdown(
    '<div class="de-section">Download Results</div>',
    unsafe_allow_html=True
)


d1, d2, d3 = st.columns(3)


with d1:

    st.download_button(
        "Download DESeq2 Results",
        data=DE_FILE.read_bytes(),
        file_name="differential_expression.csv",
        mime="text/csv",
        width="stretch"
    )


with d2:

    significant_csv = (
        significant
        .to_csv(
            index=False
        )
        .encode("utf-8")
    )

    st.download_button(
        "Download Current Significant Genes",
        data=significant_csv,
        file_name="significant_genes_current_thresholds.csv",
        mime="text/csv",
        width="stretch"
    )


with d3:

    st.download_button(
        "Download Normalized Counts",
        data=NORMALIZED_COUNTS_FILE.read_bytes(),
        file_name="normalized_counts.csv",
        mime="text/csv",
        width="stretch"
    )


# ============================================================
# FOOTER
# ============================================================

st.write("")

st.caption(
    "Method: DESeq2 • Multiple-testing correction: "
    "Benjamini-Hochberg FDR • "
    "Interactive visualizations: Plotly • "
    "PCA calculated from log2(normalized counts + 1)"
)
