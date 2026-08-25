from utils.database import get_qc_metrics

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


st.set_page_config(
    page_title="RNAFlowX | Quality Control",
    layout="wide"
)


# ============================================================
# STYLE
# ============================================================

st.markdown("""
<style>

.stApp {
    background:
        radial-gradient(circle at 10% 5%, rgba(37,99,235,.05), transparent 28%),
        radial-gradient(circle at 90% 20%, rgba(14,165,233,.04), transparent 25%),
        #F8FAFC;
}

.hero-qc {
    padding: 42px;
    border-radius: 26px;
    background: linear-gradient(135deg,#041C4A 0%,#0757C7 58%,#0EA5E9 100%);
    color: white;
    box-shadow: 0 18px 38px rgba(15,74,180,.20);
    margin-bottom: 28px;
}

.hero-qc h1 {
    font-size: 46px;
    margin: 0 0 8px 0;
    font-weight: 800;
}

.hero-qc h3 {
    color: #A5F3FC;
    margin: 0 0 15px 0;
    font-size: 21px;
    font-weight: 600;
}

.hero-qc p {
    max-width: 900px;
    color: #E8F3FF;
    font-size: 16px;
    line-height: 1.7;
}

.section-title {
    font-size: 29px;
    font-weight: 800;
    color: #123D72;
    margin: 34px 0 16px 0;
}

.metric-card {
    background: white;
    border-radius: 18px;
    padding: 21px 18px;
    border: 1px solid #E2E8F0;
    box-shadow: 0 7px 20px rgba(15,23,42,.06);
    min-height: 128px;
}

.metric-label {
    color: #64748B;
    font-size: 14px;
    font-weight: 600;
    margin-bottom: 9px;
}

.metric-value {
    color: #0F172A;
    font-size: 30px;
    font-weight: 800;
}

.metric-note {
    color: #94A3B8;
    font-size: 12px;
    margin-top: 5px;
}

.info-card {
    background: white;
    border: 1px solid #E2E8F0;
    border-radius: 18px;
    padding: 22px;
    box-shadow: 0 6px 18px rgba(15,23,42,.05);
}

</style>
""", unsafe_allow_html=True)


# ============================================================
# HERO
# ============================================================

st.markdown("""
<div class="hero-qc">

<h1>Quality Control Dashboard</h1>

<h3>
FastQC • MultiQC • Interactive Sequencing Assessment
</h3>

<p>
Interactive inspection of RNAFlowX sequencing-quality metrics derived
directly from the packaged MultiQC dataset. Metrics shown here are
generated from pipeline outputs rather than manually entered values.
</p>

</div>
""", unsafe_allow_html=True)


# ============================================================
# LOAD DATA
# ============================================================

qc = get_qc_metrics().copy()

if qc.empty:
    st.error("No MultiQC/FastQC quality-control data are available.")
    st.stop()


# ============================================================
# CLEAN DATA
# ============================================================

numeric_columns = [
    "Total Sequences",
    "Sequences flagged as poor quality",
    "Sequence length",
    "%GC",
    "total_deduplicated_percentage",
    "avg_sequence_length",
    "median_sequence_length"
]

for column in numeric_columns:
    if column in qc.columns:
        qc[column] = pd.to_numeric(
            qc[column],
            errors="coerce"
        )


def identify_stage(name):
    """
    Attempt to identify raw versus trimmed FastQC records
    using the filename/sample naming convention.
    """
    text = str(name).lower()

    trimmed_terms = [
        "trim",
        "trimmed",
        "clean",
        "fastp"
    ]

    if any(term in text for term in trimmed_terms):
        return "Trimmed"

    return "Raw"


qc["Stage"] = qc["Filename"].apply(identify_stage)

qc["Read"] = (
    qc["Sample"]
    .astype(str)
    .str.extract(r"(_[12])$", expand=False)
    .fillna("")
    .replace({
        "_1": "R1",
        "_2": "R2"
    })
)

qc["Library"] = (
    qc["Sample"]
    .astype(str)
    .str.replace(r"_[12]$", "", regex=True)
)


# ============================================================
# KPI METRICS
# ============================================================

total_fastqc_records = len(qc)

unique_libraries = qc["Library"].nunique()

total_sequences = qc["Total Sequences"].sum()

average_gc = qc["%GC"].mean()

median_read_length = qc["median_sequence_length"].median()

poor_quality = qc[
    "Sequences flagged as poor quality"
].sum()


st.markdown(
    '<div class="section-title">Sequencing Overview</div>',
    unsafe_allow_html=True
)

c1, c2, c3, c4, c5 = st.columns(5)

cards = [
    (
        "Libraries",
        f"{unique_libraries}",
        "Unique sequencing libraries"
    ),
    (
        "FastQC Records",
        f"{total_fastqc_records}",
        "Read-level QC records"
    ),
    (
        "Total Sequences",
        f"{total_sequences / 1_000_000:.1f} M",
        "Across packaged QC records"
    ),
    (
        "Average GC",
        f"{average_gc:.1f}%",
        "Mean across QC records"
    ),
    (
        "Read Length",
        f"{median_read_length:.0f} bp",
        "Median sequence length"
    )
]

for column, (label, value, note) in zip(
    [c1, c2, c3, c4, c5],
    cards
):
    with column:
        st.markdown(
            f"""
            <div class="metric-card">
                <div class="metric-label">{label}</div>
                <div class="metric-value">{value}</div>
                <div class="metric-note">{note}</div>
            </div>
            """,
            unsafe_allow_html=True
        )


# ============================================================
# SAMPLE TABLE
# ============================================================

st.markdown(
    '<div class="section-title">Sample Quality Summary</div>',
    unsafe_allow_html=True
)

summary_columns = [
    "Sample",
    "Stage",
    "Read",
    "Total Sequences",
    "Sequence length",
    "%GC",
    "total_deduplicated_percentage",
    "Sequences flagged as poor quality"
]

summary_columns = [
    col for col in summary_columns
    if col in qc.columns
]

summary = qc[summary_columns].copy()

summary = summary.rename(
    columns={
        "Sequence length": "Read Length",
        "total_deduplicated_percentage":
            "Deduplicated %",
        "Sequences flagged as poor quality":
            "Poor Quality Reads"
    }
)

st.dataframe(
    summary,
    width="stretch",
    hide_index=True,
    column_config={
        "Total Sequences": st.column_config.NumberColumn(
            format="%.0f"
        ),
        "%GC": st.column_config.NumberColumn(
            format="%.1f%%"
        ),
        "Deduplicated %": st.column_config.NumberColumn(
            format="%.1f%%"
        )
    }
)


# ============================================================
# INTERACTIVE QC CHARTS
# ============================================================

st.markdown(
    '<div class="section-title">Interactive Quality Metrics</div>',
    unsafe_allow_html=True
)

left, right = st.columns(2, gap="large")


# GC CONTENT
with left:

    st.markdown("### GC Content")

    gc_fig = px.bar(
        qc,
        x="Sample",
        y="%GC",
        color="Stage",
        hover_data=[
            "Filename",
            "Total Sequences",
            "Sequence length"
        ],
        labels={
            "%GC": "GC (%)",
            "Sample": "FASTQ"
        }
    )

    gc_fig.update_layout(
        template="plotly_white",
        height=420,
        margin=dict(
            l=20,
            r=20,
            t=35,
            b=20
        ),
        legend_title_text=""
    )

    gc_fig.update_xaxes(
        tickangle=-35
    )

    st.plotly_chart(
        gc_fig,
        width="stretch"
    )


# READ COUNTS
with right:

    st.markdown("### Sequencing Depth")

    reads_fig = px.bar(
        qc,
        x="Sample",
        y="Total Sequences",
        color="Stage",
        hover_data=[
            "%GC",
            "Sequence length",
            "Filename"
        ],
        labels={
            "Total Sequences":
                "Number of sequences",
            "Sample":
                "FASTQ"
        }
    )

    reads_fig.update_layout(
        template="plotly_white",
        height=420,
        margin=dict(
            l=20,
            r=20,
            t=35,
            b=20
        ),
        legend_title_text=""
    )

    reads_fig.update_xaxes(
        tickangle=-35
    )

    st.plotly_chart(
        reads_fig,
        width="stretch"
    )


# ============================================================
# DEDUPLICATION
# ============================================================

left, right = st.columns(2, gap="large")

with left:

    st.markdown("### Sequence Deduplication")

    dedup_fig = px.scatter(
        qc,
        x="%GC",
        y="total_deduplicated_percentage",
        color="Stage",
        symbol="Read",
        hover_name="Sample",
        hover_data=[
            "Total Sequences",
            "Sequence length"
        ],
        labels={
            "%GC": "GC (%)",
            "total_deduplicated_percentage":
                "Deduplicated sequences (%)"
        }
    )

    dedup_fig.update_traces(
        marker=dict(size=12)
    )

    dedup_fig.update_layout(
        template="plotly_white",
        height=410,
        margin=dict(
            l=20,
            r=20,
            t=35,
            b=20
        )
    )

    st.plotly_chart(
        dedup_fig,
        width="stretch"
    )


# ============================================================
# POOR QUALITY
# ============================================================

with right:

    st.markdown("### Reads Flagged as Poor Quality")

    # --------------------------------------------------------
    # Calculate QC statistics
    # --------------------------------------------------------

    total_flagged = int(
        qc["Sequences flagged as poor quality"]
        .fillna(0)
        .sum()
    )

    total_reads_qc = int(
        qc["Total Sequences"]
        .fillna(0)
        .sum()
    )

    if total_reads_qc > 0:
        poor_quality_rate = (
            total_flagged / total_reads_qc
        ) * 100
    else:
        poor_quality_rate = 0.0

    # --------------------------------------------------------
    # PASS STATE
    # --------------------------------------------------------

    if total_flagged == 0:

        st.success("✓ QC PASS — No Poor-Quality Reads Detected")

        q1, q2 = st.columns(2)

        with q1:
            st.metric(
                label="Reads Flagged",
                value="0"
            )

        with q2:
            st.metric(
                label="Poor-quality Rate",
                value="0.00%"
            )

        st.markdown(
            """
**Sequencing quality assessment**

No sequences were flagged as poor quality across the analyzed
FastQC records.

This indicates that the sequencing libraries passed this
quality-control criterion without poor-quality read flags.
"""
        )

        st.progress(100)

        st.caption(
            f"{total_reads_qc:,} sequences evaluated • "
            "FastQC poor-quality sequence flag: PASS"
        )

    # --------------------------------------------------------
    # FLAGGED READS STATE
    # --------------------------------------------------------

    else:

        st.warning(
            f"{total_flagged:,} sequences were flagged "
            "as poor quality."
        )

        q1, q2 = st.columns(2)

        with q1:
            st.metric(
                label="Reads Flagged",
                value=f"{total_flagged:,}"
            )

        with q2:
            st.metric(
                label="Poor-quality Rate",
                value=f"{poor_quality_rate:.3f}%"
            )

        poor_plot_df = qc[
            [
                "Sample",
                "Sequences flagged as poor quality"
            ]
        ].copy()

        poor_plot_df[
            "Sequences flagged as poor quality"
        ] = pd.to_numeric(
            poor_plot_df[
                "Sequences flagged as poor quality"
            ],
            errors="coerce"
        ).fillna(0)

        poor_fig = px.scatter(
            poor_plot_df,
            x="Sample",
            y="Sequences flagged as poor quality",
            size="Sequences flagged as poor quality",
            color="Sequences flagged as poor quality",
            color_continuous_scale="Reds",
            hover_name="Sample",
            labels={
                "Sample": "FASTQ",
                "Sequences flagged as poor quality":
                    "Flagged sequences"
            }
        )

        poor_fig.update_traces(
            marker=dict(
                sizemin=8,
                line=dict(
                    width=1
                )
            )
        )

        poor_fig.update_layout(
            template="plotly_white",
            height=340,
            coloraxis_colorbar=dict(
                title="Flagged"
            ),
            margin=dict(
                l=20,
                r=20,
                t=20,
                b=20
            )
        )

        st.plotly_chart(
            poor_fig,
            width="stretch"
        )

        st.caption(
            f"{total_reads_qc:,} sequences evaluated • "
            f"{poor_quality_rate:.3f}% flagged"
        )

# ============================================================
# FASTQC MODULE STATUS
# ============================================================

st.markdown(
    '<div class="section-title">FastQC Module Performance</div>',
    unsafe_allow_html=True
)

module_columns = [
    "basic_statistics",
    "per_base_sequence_quality",
    "per_tile_sequence_quality",
    "per_sequence_quality_scores",
    "per_base_sequence_content",
    "per_sequence_gc_content",
    "per_base_n_content",
    "sequence_length_distribution",
    "sequence_duplication_levels",
    "overrepresented_sequences",
    "adapter_content"
]

module_columns = [
    col for col in module_columns
    if col in qc.columns
]

status_records = []

for module in module_columns:

    counts = (
        qc[module]
        .fillna("unknown")
        .astype(str)
        .str.lower()
        .value_counts()
    )

    for status in [
        "pass",
        "warn",
        "fail",
        "unknown"
    ]:

        count = counts.get(status, 0)

        if count > 0:
            status_records.append({
                "Module":
                    module.replace("_", " ").title(),
                "Status":
                    status.upper(),
                "Count":
                    count
            })

status_df = pd.DataFrame(status_records)

if not status_df.empty:

    module_fig = px.bar(
        status_df,
        x="Module",
        y="Count",
        color="Status",
        barmode="stack",
        category_orders={
            "Status": [
                "PASS",
                "WARN",
                "FAIL",
                "UNKNOWN"
            ]
        }
    )

    module_fig.update_layout(
        template="plotly_white",
        height=470,
        margin=dict(
            l=20,
            r=20,
            t=30,
            b=100
        ),
        legend_title_text="FastQC status"
    )

    module_fig.update_xaxes(
        tickangle=-35
    )

    st.plotly_chart(
        module_fig,
        width="stretch"
    )


# ============================================================
# QC STATUS TABLE
# ============================================================

status_matrix = qc[
    ["Sample"] + module_columns
].copy()

status_matrix.columns = [
    column.replace("_", " ").title()
    for column in status_matrix.columns
]

st.dataframe(
    status_matrix,
    width="stretch",
    hide_index=True
)


# ============================================================
# INTERPRETATION
# ============================================================

st.markdown(
    '<div class="section-title">Quality-Control Interpretation</div>',
    unsafe_allow_html=True
)

left, right = st.columns(2, gap="large")

with left:

    pass_count = sum(
        (qc[col].astype(str).str.lower() == "pass").sum()
        for col in module_columns
    )

    warning_count = sum(
        (qc[col].astype(str).str.lower() == "warn").sum()
        for col in module_columns
    )

    fail_count = sum(
        (qc[col].astype(str).str.lower() == "fail").sum()
        for col in module_columns
    )

    st.markdown(
        f"""
        <div class="info-card">
        <h3>QC Status</h3>

        <p><strong>PASS:</strong> {pass_count}</p>
        <p><strong>WARN:</strong> {warning_count}</p>
        <p><strong>FAIL:</strong> {fail_count}</p>

        <p>
        FastQC warnings and failures are displayed rather than hidden,
        allowing sequencing-quality issues to remain transparent during
        downstream interpretation.
        </p>

        </div>
        """,
        unsafe_allow_html=True
    )

with right:

    st.markdown(
        f"""
        <div class="info-card">
        <h3>Dataset Characteristics</h3>

        <p><strong>Average GC:</strong> {average_gc:.1f}%</p>
        <p><strong>Median read length:</strong> {median_read_length:.0f} bp</p>
        <p><strong>Poor-quality reads flagged:</strong> {poor_quality:,.0f}</p>
        <p><strong>FastQC records evaluated:</strong> {total_fastqc_records}</p>

        <p>
        Interactive plots support hover inspection, zooming,
        filtering through Plotly controls and sample-level comparison.
        </p>

        </div>
        """,
        unsafe_allow_html=True
    )


# ============================================================
# SOURCE
# ============================================================

st.write("")

st.caption(
    "Source: RNAFlowX packaged MultiQC/FastQC database. "
    "Displayed values are derived from pipeline-generated QC outputs."
)
