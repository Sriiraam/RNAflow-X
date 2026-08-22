import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="RNAFlowX | Downloads",
    layout="wide"
)

# ============================================================
# PATHS
# ============================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS = PROJECT_ROOT / "results"

# ============================================================
# CSS
# ============================================================
st.markdown("""
<style>

.download-hero{
    padding:46px;
    border-radius:26px;
    background:linear-gradient(135deg,#062A67 0%,#0E56D8 100%);
    color:white;
    box-shadow:0 18px 36px rgba(15,74,180,.20);
    margin-bottom:30px;
}

.download-hero h1{
    margin:0 0 12px;
    font-size:48px;
    font-weight:800;
}

.download-hero h3{
    margin:0 0 18px;
    color:#75E6FF;
    font-size:22px;
}

.download-hero p{
    max-width:900px;
    color:#E8F3FF;
    font-size:16px;
    line-height:1.8;
}

.section-title{
    margin:34px 0 16px;
    font-size:30px;
    font-weight:800;
    color:#153F73;
}

.download-card{
    background:white;
    border-radius:18px;
    padding:22px;
    box-shadow:0 7px 20px rgba(15,23,42,.07);
    border-top:5px solid #2563EB;
    min-height:170px;
}

.download-card h3{
    color:#153F73;
    margin-bottom:10px;
}

.download-card p{
    color:#64748B;
    line-height:1.6;
}

.status-ok{
    color:#15803D;
    font-weight:700;
}

.status-missing{
    color:#DC2626;
    font-weight:700;
}

</style>
""", unsafe_allow_html=True)

# ============================================================
# HELPERS
# ============================================================
def render_download_card(title, description, file_path, mime, key):
    st.markdown(f"""
<div class="download-card">
<h3>{title}</h3>
<p>{description}</p>
</div>
""", unsafe_allow_html=True)

    if file_path.exists():
        st.markdown(
            '<div class="status-ok">✓ Available</div>',
            unsafe_allow_html=True
        )

        with open(file_path, "rb") as f:
            st.download_button(
                label=f"Download {title}",
                data=f,
                file_name=file_path.name,
                mime=mime,
                width="stretch",
                key=key
            )
    else:
        st.markdown(
            '<div class="status-missing">✗ File not found</div>',
            unsafe_allow_html=True
        )


# ============================================================
# HERO
# ============================================================
st.markdown("""
<div class="download-hero">

<h1>Downloads & Scientific Reports</h1>

<h3>
Analysis Outputs • Tables • Reports • Reproducible Results
</h3>

<p>
Download the principal RNAFlowX outputs generated across quantification,
differential expression, functional enrichment and quality-control stages.
These files provide reusable inputs for downstream analysis, reporting and
independent validation.
</p>

</div>
""", unsafe_allow_html=True)

# ============================================================
# DIFFERENTIAL EXPRESSION
# ============================================================
st.markdown(
    '<div class="section-title">🧬 Differential Expression</div>',
    unsafe_allow_html=True
)

de_dir = RESULTS / "differential_expression" / "deseq2_results"

c1, c2, c3 = st.columns(3, gap="large")

with c1:
    render_download_card(
        "Differential Expression",
        "Complete DESeq2 results including baseMean, log2 fold change, p-value and adjusted p-value.",
        de_dir / "differential_expression.csv",
        "text/csv",
        "de_all"
    )

with c2:
    render_download_card(
        "Significant Genes",
        "Filtered significant differential-expression results for downstream interpretation.",
        de_dir / "significant_genes.csv",
        "text/csv",
        "de_sig"
    )

with c3:
    render_download_card(
        "Normalized Counts",
        "DESeq2-normalized expression matrix across all biological samples.",
        de_dir / "normalized_counts.csv",
        "text/csv",
        "de_norm"
    )

# ============================================================
# COUNTING
# ============================================================
st.markdown(
    '<div class="section-title">📊 Quantification</div>',
    unsafe_allow_html=True
)

q1, q2 = st.columns(2, gap="large")

with q1:
    render_download_card(
        "Count Matrix",
        "Gene-level count matrix used as input for differential-expression analysis.",
        RESULTS / "counting" / "count_matrix.csv",
        "text/csv",
        "count_matrix"
    )

with q2:
    render_download_card(
        "TxImport Summary",
        "Summary generated during transcript-to-gene quantification import.",
        RESULTS / "counting" / "tximport_summary.txt",
        "text/plain",
        "tximport"
    )

# ============================================================
# FUNCTIONAL ENRICHMENT
# ============================================================
st.markdown(
    '<div class="section-title">🔬 Functional Enrichment</div>',
    unsafe_allow_html=True
)

e1, e2, e3 = st.columns(3, gap="large")

with e1:
    render_download_card(
        "GO Biological Process",
        "Gene Ontology Biological Process enrichment for all significant genes.",
        RESULTS / "enrichment" / "GO" / "BP" / "all_significant.csv",
        "text/csv",
        "go_bp"
    )

with e2:
    render_download_card(
        "KEGG Pathways",
        "KEGG pathway enrichment for all significant genes.",
        RESULTS / "enrichment" / "KEGG" / "all_significant.csv",
        "text/csv",
        "kegg"
    )

with e3:
    render_download_card(
        "GSEA Results",
        "GO Biological Process Gene Set Enrichment Analysis results.",
        RESULTS / "enrichment" / "GSEA" / "GO_BP_GSEA.csv",
        "text/csv",
        "gsea"
    )

# ============================================================
# SUMMARIES
# ============================================================
st.markdown(
    '<div class="section-title">📄 Analysis Summaries</div>',
    unsafe_allow_html=True
)

s1, s2 = st.columns(2, gap="large")

with s1:
    render_download_card(
        "DESeq2 Summary",
        "Text summary of the differential-expression analysis.",
        de_dir / "deseq2_summary.txt",
        "text/plain",
        "deseq_summary"
    )

with s2:
    render_download_card(
        "Enrichment Summary",
        "Summary of GO, KEGG and GSEA enrichment analyses.",
        RESULTS / "enrichment" / "enrichment_summary.txt",
        "text/plain",
        "enrichment_summary"
    )

# ============================================================
# MULTIQC
# ============================================================
st.markdown(
    '<div class="section-title">🧪 Quality Control Report</div>',
    unsafe_allow_html=True
)

render_download_card(
    "MultiQC Report",
    "Complete interactive MultiQC HTML report summarizing sequencing-quality metrics across samples.",
    RESULTS / "multiqc" / "multiqc_report.html",
    "text/html",
    "multiqc_report"
)

# ============================================================
# R OBJECT
# ============================================================
st.markdown(
    '<div class="section-title">⚙️ Reproducibility Artifact</div>',
    unsafe_allow_html=True
)

render_download_card(
    "DESeq2 Dataset Object",
    "Serialized DESeq2 dataset containing the fitted analysis object for reproducibility and downstream R-based investigation.",
    de_dir / "dds.rds",
    "application/octet-stream",
    "dds_rds"
)

st.write("")
st.caption(
    "RNAFlowX downloads are generated directly from pipeline outputs. "
    "No example or placeholder files are used."
)
