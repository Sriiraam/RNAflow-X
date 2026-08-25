from utils.database import get_go, get_kegg

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go


# ============================================================
# PAGE CONFIG
# ============================================================

st.set_page_config(
    page_title="RNAFlowX | Functional Enrichment",
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
        rgba(37,99,235,.045), transparent 28%),
        #F8FAFC;
}

.enrich-hero{
    padding:42px 44px;
    border-radius:24px;
    background:
        linear-gradient(
            135deg,
            #031A46 0%,
            #064D9E 55%,
            #0877C9 100%
        );
    color:white;
    box-shadow:0 18px 38px rgba(15,74,180,.18);
    margin-bottom:24px;
}

.enrich-hero h1{
    margin:0;
    font-size:46px;
    font-weight:850;
}

.enrich-hero h3{
    margin:8px 0 14px;
    color:#5EEAD4;
    font-size:20px;
}

.enrich-hero p{
    max-width:920px;
    color:#E8F3FF;
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
    background:#EFF6FF;
    border-left:5px solid #2563EB;
    padding:20px 22px;
    border-radius:14px;
    color:#334155;
    line-height:1.75;
}

</style>
""", unsafe_allow_html=True)


# ============================================================
# HELPERS
# ============================================================

def ratio_to_float(value):

    try:
        numerator, denominator = str(value).split("/")
        return float(numerator) / float(denominator)

    except Exception:
        return np.nan


def prepare_enrichment(df):

    df = df.copy()

    if df.empty:
        return df

    for col in [
        "pvalue",
        "p.adjust",
        "qvalue",
        "Count"
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col],
                errors="coerce"
            )

    if "GeneRatio" in df.columns:
        df["GeneRatio_numeric"] = (
            df["GeneRatio"]
            .apply(ratio_to_float)
        )

    df["p_adjust_plot"] = (
        df["p.adjust"]
        .fillna(1)
        .clip(lower=1e-300)
    )

    df["neg_log10_padj"] = (
        -np.log10(
            df["p_adjust_plot"]
        )
    )

    return df


def top_terms(df, n=15):

    return (
        df
        .sort_values(
            "p.adjust",
            ascending=True
        )
        .head(n)
        .copy()
    )


# ============================================================
# LOAD DEFAULT DATA
# ============================================================

go_bp_all = prepare_enrichment(
    get_go("BP", "all")
)

kegg_all = prepare_enrichment(
    get_kegg("all")
)


# ============================================================
# HERO
# ============================================================

st.markdown("""
<div class="enrich-hero">

<h1>Functional Enrichment Analysis</h1>

<h3>
Gene Ontology • KEGG • Interactive Pathway Exploration
</h3>

<p>
Functional enrichment converts differentially expressed genes into
biologically interpretable processes and pathways. RNAFlowX integrates
Gene Ontology and KEGG enrichment to identify coordinated biological
responses associated with PFOS exposure.
</p>

</div>
""", unsafe_allow_html=True)


# ============================================================
# KPI
# ============================================================

go_sig = go_bp_all[
    go_bp_all["p.adjust"] < 0.05
]

kegg_sig = kegg_all[
    kegg_all["p.adjust"] < 0.05
]

k1, k2, k3, k4 = st.columns(4)

cards = [
    (
        k1,
        "GO BP Terms",
        len(go_sig),
        "#2563EB",
        "Adjusted p < 0.05"
    ),
    (
        k2,
        "KEGG Pathways",
        len(kegg_sig),
        "#16A34A",
        "Adjusted p < 0.05"
    ),
    (
        k3,
        "Top GO Signal",
        f"{go_sig['neg_log10_padj'].max():.1f}"
        if not go_sig.empty else "0",
        "#7C3AED",
        "-log10 adjusted p"
    ),
    (
        k4,
        "Top KEGG Signal",
        f"{kegg_sig['neg_log10_padj'].max():.1f}"
        if not kegg_sig.empty else "0",
        "#EA580C",
        "-log10 adjusted p"
    )
]

for col, label, value, color, note in cards:

    with col:

        st.markdown(
            f"""
<div class="kpi" style="border-top:5px solid {color};">
<div class="kpi-label">{label}</div>
<div class="kpi-value" style="color:{color};">{value}</div>
<div class="kpi-note">{note}</div>
</div>
""",
            unsafe_allow_html=True
        )


# ============================================================
# EXPLORER
# ============================================================

st.markdown(
    '<div class="section-title">Enrichment Explorer</div>',
    unsafe_allow_html=True
)

tab_go, tab_kegg = st.tabs(
    [
        "Gene Ontology",
        "KEGG Pathways"
    ]
)


# ============================================================
# GO
# ============================================================

with tab_go:

    c1, c2, c3 = st.columns(
        [1.2, 1.2, 1]
    )

    with c1:

        go_type = st.selectbox(
            "Ontology",
            [
                "Biological Process",
                "Molecular Function",
                "Cellular Component"
            ]
        )

    with c2:

        regulation = st.selectbox(
            "Gene set",
            [
                "All Significant",
                "Upregulated",
                "Downregulated"
            ],
            key="go_regulation"
        )

    with c3:

        top_n = st.slider(
            "Terms displayed",
            min_value=5,
            max_value=30,
            value=15,
            step=5
        )


    go_code = {
        "Biological Process": "BP",
        "Molecular Function": "MF",
        "Cellular Component": "CC"
    }[go_type]

    regulation_code = {
        "All Significant": "all",
        "Upregulated": "up",
        "Downregulated": "down"
    }[regulation]


    go_df = prepare_enrichment(
        get_go(
            go_code,
            regulation_code
        )
    )


    if go_df.empty:

        st.warning(
            "No GO enrichment results "
            "are available for this selection."
        )

    else:

        go_sig_selected = go_df[
            go_df["p.adjust"] < 0.05
        ].copy()

        if go_sig_selected.empty:

            st.info(
                "No GO terms meet adjusted p < 0.05."
            )

        else:

            top_go = top_terms(
                go_sig_selected,
                top_n
            )


            # ================================================
            # DOT + BAR
            # ================================================

            st.markdown(
                '<div class="section-title">'
                'GO Enrichment Landscape'
                '</div>',
                unsafe_allow_html=True
            )

            left, right = st.columns(
                2,
                gap="large"
            )


            with left:

                st.markdown(
                    "### Interactive GO Dot Plot"
                )

                go_dot = px.scatter(
                    top_go.sort_values(
                        "GeneRatio_numeric"
                    ),

                    x="GeneRatio_numeric",
                    y="Description",

                    size="Count",
                    color="neg_log10_padj",

                    color_continuous_scale="Turbo",

                    hover_data={
                        "GeneRatio": True,
                        "Count": True,
                        "p.adjust": ":.2e",
                        "pvalue": ":.2e",
                        "GeneRatio_numeric": False,
                        "neg_log10_padj": ":.2f"
                    },

                    labels={
                        "GeneRatio_numeric":
                            "Gene Ratio",

                        "neg_log10_padj":
                            "-log10(padj)",

                        "Description":
                            ""
                    }
                )

                go_dot.update_layout(
                    template="plotly_white",
                    height=560,
                    margin=dict(
                        l=10,
                        r=10,
                        t=20,
                        b=20
                    )
                )

                st.plotly_chart(
                    go_dot,
                    width="stretch"
                )


            with right:

                st.markdown(
                    "### Ranked Biological Processes"
                )

                go_bar = px.bar(
                    top_go.sort_values(
                        "neg_log10_padj"
                    ),

                    x="neg_log10_padj",
                    y="Description",

                    orientation="h",

                    color="GeneRatio_numeric",

                    color_continuous_scale="Viridis",

                    hover_data={
                        "Count": True,
                        "GeneRatio": True,
                        "p.adjust": ":.2e"
                    },

                    labels={
                        "neg_log10_padj":
                            "-log10 adjusted p-value",

                        "GeneRatio_numeric":
                            "Gene Ratio",

                        "Description":
                            ""
                    }
                )

                go_bar.update_layout(
                    template="plotly_white",
                    height=560,
                    margin=dict(
                        l=10,
                        r=10,
                        t=20,
                        b=20
                    )
                )

                st.plotly_chart(
                    go_bar,
                    width="stretch"
                )


            # ================================================
            # JOINT STYLE
            # ================================================

            st.markdown(
                '<div class="section-title">'
                'Gene Ratio & Significance Relationship'
                '</div>',
                unsafe_allow_html=True
            )

            joint_fig = px.scatter(
            go_sig_selected,

            x="GeneRatio_numeric",
            y="neg_log10_padj",

            size="Count",

            color="neg_log10_padj",

            color_continuous_scale="Plasma",

            hover_name="Description",

            hover_data={
            "GeneRatio": True,
            "Count": True,
            "p.adjust": ":.2e",
            "GeneRatio_numeric": False,
            "neg_log10_padj": ":.2f"
            },

            labels={
            "GeneRatio_numeric":
            "Gene Ratio",

            "neg_log10_padj":
            "-log10 adjusted p-value"
            }
            )

            joint_fig.update_traces(
            marker=dict(
            opacity=0.78,
            line=dict(
            width=0.7,
            color="white"
            )
            )
            )

            joint_fig.update_layout(
            template="plotly_white",
            height=540,
            margin=dict(
            l=25,
            r=20,
            t=20,
            b=25
            ),
            coloraxis_colorbar=dict(
            title="-log10(padj)"
            )
            )

            st.plotly_chart(
            joint_fig,
            width="stretch"
            )


            # ================================================
            # TABLE
            # ================================================

            st.markdown(
                '<div class="section-title">'
                'GO Term Explorer'
                '</div>',
                unsafe_allow_html=True
            )

            search = st.text_input(
                "Search GO term",
                placeholder="Example: blood vessel development"
            )

            display_go = (
                go_sig_selected.copy()
            )

            if search:

                display_go = display_go[
                    display_go["Description"]
                    .astype(str)
                    .str.contains(
                        search,
                        case=False,
                        na=False
                    )
                ]

            display_go = (
                display_go
                .sort_values(
                    "p.adjust"
                )
            )

            columns = [
                "ID",
                "Description",
                "GeneRatio",
                "BgRatio",
                "Count",
                "pvalue",
                "p.adjust",
                "qvalue"
            ]

            st.dataframe(
                display_go[
                    [
                        c for c in columns
                        if c in display_go.columns
                    ]
                ],
                width="stretch",
                hide_index=True,
                height=430,

                column_config={
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


# ============================================================
# KEGG
# ============================================================

with tab_kegg:

    c1, c2 = st.columns(
        [1.5, 1]
    )

    with c1:

        kegg_regulation = st.selectbox(
            "Gene set",
            [
                "All Significant",
                "Upregulated",
                "Downregulated"
            ],
            key="kegg_reg"
        )

    with c2:

        kegg_top_n = st.slider(
            "Pathways displayed",
            5,
            20,
            10,
            key="kegg_top"
        )


    regulation_code = {
        "All Significant": "all",
        "Upregulated": "up",
        "Downregulated": "down"
    }[kegg_regulation]


    kegg_df = prepare_enrichment(
        get_kegg(
            regulation_code
        )
    )


    if kegg_df.empty:

        st.warning(
            "No KEGG results available."
        )

    else:

        kegg_sig_selected = kegg_df[
            kegg_df["p.adjust"] < 0.05
        ].copy()

        top_kegg = top_terms(
            kegg_sig_selected,
            kegg_top_n
        )


        if top_kegg.empty:

            st.info(
                "No KEGG pathways meet "
                "adjusted p < 0.05."
            )

        else:

            st.markdown(
                '<div class="section-title">'
                'KEGG Pathway Landscape'
                '</div>',
                unsafe_allow_html=True
            )

            left, right = st.columns(
                2,
                gap="large"
            )


            with left:

                kegg_dot = px.scatter(
                    top_kegg.sort_values(
                        "GeneRatio_numeric"
                    ),

                    x="GeneRatio_numeric",
                    y="Description",

                    size="Count",
                    color="neg_log10_padj",

                    color_continuous_scale="Turbo",

                    hover_data={
                        "GeneRatio": True,
                        "Count": True,
                        "p.adjust": ":.2e",
                        "GeneRatio_numeric": False
                    },

                    labels={
                        "GeneRatio_numeric":
                            "Gene Ratio",

                        "neg_log10_padj":
                            "-log10(padj)",

                        "Description":
                            ""
                    }
                )

                kegg_dot.update_layout(
                    template="plotly_white",
                    height=520
                )

                st.plotly_chart(
                    kegg_dot,
                    width="stretch"
                )


            with right:

                kegg_bar = px.bar(
                    top_kegg.sort_values(
                        "neg_log10_padj"
                    ),

                    x="neg_log10_padj",
                    y="Description",

                    orientation="h",

                    color="Count",

                    color_continuous_scale="Sunset",

                    hover_data={
                        "GeneRatio": True,
                        "p.adjust": ":.2e"
                    },

                    labels={
                        "neg_log10_padj":
                            "-log10 adjusted p-value",

                        "Description":
                            ""
                    }
                )

                kegg_bar.update_layout(
                    template="plotly_white",
                    height=520
                )

                st.plotly_chart(
                    kegg_bar,
                    width="stretch"
                )


            # Histogram
            st.markdown(
                '<div class="section-title">'
                'KEGG Significance Distribution'
                '</div>',
                unsafe_allow_html=True
            )

            hist_fig = px.histogram(
                kegg_sig_selected,

                x="neg_log10_padj",

                nbins=12,

                color_discrete_sequence=[
                    "#7C3AED"
                ],

                hover_data=[
                    "Description"
                ],

                labels={
                    "neg_log10_padj":
                        "-log10 adjusted p-value"
                }
            )

            hist_fig.update_layout(
                template="plotly_white",
                height=390,
                yaxis_title="Number of pathways"
            )

            st.plotly_chart(
                hist_fig,
                width="stretch"
            )


            # Table
            st.markdown(
                '<div class="section-title">'
                'KEGG Pathway Explorer'
                '</div>',
                unsafe_allow_html=True
            )

            search_kegg = st.text_input(
                "Search KEGG pathway",
                placeholder="Example: focal adhesion"
            )

            display_kegg = (
                kegg_sig_selected.copy()
            )

            if search_kegg:

                display_kegg = display_kegg[
                    display_kegg["Description"]
                    .astype(str)
                    .str.contains(
                        search_kegg,
                        case=False,
                        na=False
                    )
                ]

            columns = [
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

            st.dataframe(
                display_kegg[
                    [
                        c for c in columns
                        if c in display_kegg.columns
                    ]
                ].sort_values(
                    "p.adjust"
                ),

                width="stretch",
                hide_index=True,
                height=430,

                column_config={
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


# ============================================================
# INTERPRETATION
# ============================================================

st.markdown(
    '<div class="section-title">'
    'Biological Interpretation'
    '</div>',
    unsafe_allow_html=True
)

if not go_sig.empty:

    strongest_go = (
        go_sig
        .sort_values("p.adjust")
        .iloc[0]
    )

    go_text = strongest_go[
        "Description"
    ]

else:
    go_text = "No significant GO process"


if not kegg_sig.empty:

    strongest_kegg = (
        kegg_sig
        .sort_values("p.adjust")
        .iloc[0]
    )

    kegg_text = strongest_kegg[
        "Description"
    ]

else:
    kegg_text = "No significant KEGG pathway"


st.markdown(
    f"""
<div class="info-box">

<strong>Strongest GO Biological Process:</strong>
{go_text}

<br><br>

<strong>Strongest KEGG pathway:</strong>
{kegg_text}

<br><br>

Gene Ontology identifies enriched biological functions,
whereas KEGG places differentially expressed genes into
curated molecular and signaling pathways. GeneRatio,
gene count and multiple-testing-adjusted significance
should be interpreted together rather than relying on
p-values alone.

</div>
""",
    unsafe_allow_html=True
)


st.write("")

st.caption(
    "Functional enrichment source: clusterProfiler • "
    "Significance criterion: adjusted p-value < 0.05"
)
