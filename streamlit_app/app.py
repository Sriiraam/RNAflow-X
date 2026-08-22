import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="RNAFlowX",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ---------- Load CSS ----------

st.markdown("""
<style>
.main > div{
    background:#F4F8FC;
}

/* ================= HERO ================= */
.hero{
position:relative;
padding:52px 46px 60px;
margin-bottom:34px;
border-radius:28px;
overflow:hidden;
background:linear-gradient(135deg,#062A67 0%,#0E56D8 100%);
box-shadow:0 20px 38px rgba(16,74,180,.20);
color:white;
min-height:315px;
}

.hero:before{
content:"";
position:absolute;
right:-120px;
bottom:-100px;
width:480px;
height:240px;
background:radial-gradient(circle,rgba(255,255,255,.12),transparent 72%);
}

.hero h1{
font-size:60px;
margin:18px 0 12px;
font-weight:800;
letter-spacing:-1px;
}

.hero h3{
font-size:26px;
margin-bottom:22px;
font-weight:600;
color:#73E7FF;
}

.hero p{
font-size:18px;
line-height:1.9;
max-width:900px;
color:#E8F3FF;
}

.hero small{
font-size:13px;
letter-spacing:2px;
font-weight:700;
color:#B9D8FF;
}

/* glowing dots */
.dot{
position:absolute;
border-radius:50%;
background:#66E5FF;
box-shadow:0 0 18px #66E5FF;
opacity:.9;
}
.d1{width:14px;height:14px;top:40px;right:130px;}
.d2{width:8px;height:8px;top:120px;right:240px;}
.d3{width:18px;height:18px;top:165px;right:70px;}

/* ================= KPI ================= */
.kpi{
background:white;
border-radius:18px;
padding:22px 18px;
text-align:center;
box-shadow:0 8px 20px rgba(15,23,42,.08);
border-top:6px solid;
height:150px;
display:flex;
flex-direction:column;
justify-content:center;
}

.kpi .icon{
font-size:30px;
margin-bottom:8px;
}

.kpi .value{
font-size:42px;
font-weight:800;
line-height:1;
margin-bottom:10px;
color:#0F172A;
}

.kpi .label{
font-size:15px;
font-weight:600;
color:#64748B;
}

/* ================= CARDS ================= */
.card{
background:white;
padding:30px;
border-radius:22px;
box-shadow:0 8px 20px rgba(15,23,42,.06);
height:100%;
}

.card h2{
margin-bottom:22px;
font-size:32px;
color:#144B8B;
}

.card p{
font-size:17px;
line-height:1.9;
color:#475569;
}

/* output box */
.outputs{
background:#EFF6FF;
border-left:5px solid #2563EB;
padding:16px;
border-radius:12px;
margin-top:18px;
}

/* ================= MINI ================= */
.mini{
background:#F8FAFC;
border-left:5px solid;
border-radius:14px;
padding:14px;
margin-bottom:12px;
}

.mini .t{
font-size:11px;
font-weight:700;
letter-spacing:.8px;
color:#64748B;
text-transform:uppercase;
}

.mini .v{
font-size:24px;
font-weight:800;
margin-top:6px;
color:#0F172A;
}

.card-center{
background:white;
border-radius:18px;
padding:22px 18px;
text-align:center;
box-shadow:0 6px 18px rgba(15,23,42,.06);
min-height:210px;
display:flex;
flex-direction:column;
justify-content:flex-start;
align-items:center;
}

.card-center h4{
margin:12px 0 8px;
font-size:22px;
font-weight:700;
color:#144B8B;
}

.card-center p{
margin:0;
font-size:14px;
line-height:1.6;
color:#64748B;
}
</style>
""", unsafe_allow_html=True)

# ---------- Sidebar ----------

st.sidebar.markdown("# 🧬 RNAFlowX")
st.sidebar.caption("Bulk RNA-seq Engineering Platform")

st.sidebar.success("✅ Pipeline Completed")

st.sidebar.divider()

st.sidebar.markdown("### 🛠 Tech Stack")
st.sidebar.markdown("""
- Nextflow DSL2
- FastQC / MultiQC
- STAR Aligner
- DESeq2
- clusterProfiler
- SQLite
- Streamlit
- Docker
""")

# ---------- HERO ----------

st.markdown("""
<div class="hero">

<div class="dot d1"></div>
<div class="dot d2"></div>
<div class="dot d3"></div>

<small>TRANSCRIPTOMICS • WORKFLOW ENGINEERING • REPRODUCIBLE BIOINFORMATICS</small>

<h1>RNAFlowX</h1>

<h3>Production-grade Bulk RNA Sequencing Analysis Platform</h3>

<p>
RNAFlowX transforms raw paired-end FASTQ sequencing reads into biologically
interpretable transcriptomic insights using a modular Nextflow DSL2 workflow.
The platform integrates quality assessment, transcript quantification,
differential expression, pathway enrichment and interactive visualization while
maintaining complete computational reproducibility.
</p>

</div>
""", unsafe_allow_html=True)

st.write("")

# ---------- KPI ----------

k1,k2,k3,k4 = st.columns(4,gap="medium")

with k1:
    st.markdown("""
<div class="kpi" style="border-color:#2563EB">
<div class="icon">🧬</div>
<div class="value">10,283</div>
<div class="label">Genes Analysed</div>
</div>
""",unsafe_allow_html=True)

with k2:
    st.markdown("""
<div class="kpi" style="border-color:#16A34A">
<div class="icon">📈</div>
<div class="value">430</div>
<div class="label">Significant DEGs</div>
</div>
""",unsafe_allow_html=True)

with k3:
    st.markdown("""
<div class="kpi" style="border-color:#7C3AED">
<div class="icon">👥</div>
<div class="value">4</div>
<div class="label">Biological Samples</div>
</div>
""",unsafe_allow_html=True)

with k4:
    st.markdown("""
<div class="kpi" style="border-color:#EA580C">
<div class="icon">⚙️</div>
<div class="value">9</div>
<div class="label">Pipeline Modules</div>
</div>
""",unsafe_allow_html=True)

st.write("")

# ---------- Executive Overview ----------

left,right=st.columns([1.35,1],gap="large")

with left:
    st.markdown("""
<div class="card">

<h2>📋 Executive Project Overview</h2>

<p>
RNAFlowX is a modular transcriptomics platform engineered to demonstrate
modern bioinformatics workflow development rather than a single RNA-seq analysis.
</p>

<p>
The workflow separates every analytical stage into reusable DSL2 modules,
enabling identical execution across local workstations, Docker containers,
HPC clusters and future cloud environments through configurable execution profiles.
</p>

<div class="outputs">
<b>Primary Outputs</b><br>
FastQC • STAR Alignment • featureCounts • DESeq2 • GO • KEGG • GSEA • MultiQC
</div>

</div>
""",unsafe_allow_html=True)

with right:
    st.markdown("""
<div class="card">

<h2>🧪 Experimental Design</h2>

<div class="mini" style="border-color:#2563EB">
<div class="t">Organism</div>
<div class="v">H. sapiens</div>
</div>

<div class="mini" style="border-color:#16A34A">
<div class="t">Condition</div>
<div class="v">PFOS</div>
</div>

<div class="mini" style="border-color:#7C3AED">
<div class="t">Controls</div>
<div class="v">2</div>
</div>

<div class="mini" style="border-color:#0EA5E9">
<div class="t">Treated</div>
<div class="v">2</div>
</div>

<div class="mini" style="border-color:#EA580C">
<div class="t">Read Type</div>
<div class="v">Paired-end</div>
</div>

<div class="mini" style="border-color:#DC2626">
<div class="t">Read Length</div>
<div class="v">36 bp</div>
</div>

</div>
""",unsafe_allow_html=True)

st.write("")

# ---------- Pipeline Architecture ----------
st.markdown("## RNAFlowX Pipeline Architecture")

p1, p2, p3, p4, p5 = st.columns(5)

with p1:
    st.markdown("""
<div class="card-center">
<div style="font-size:32px;">01</div>
<h4>FASTQ</h4>
<p>Raw sequencing reads</p>
</div>
""", unsafe_allow_html=True)

with p2:
    st.markdown("""
<div class="card-center">
<div style="font-size:32px;">02</div>
<h4>QC & Trim</h4>
<p>FastQC + FastP</p>
</div>
""", unsafe_allow_html=True)

with p3:
    st.markdown("""
<div class="card-center">
<div style="font-size:32px;">03</div>
<h4>Alignment</h4>
<p>STAR + SAMtools</p>
</div>
""", unsafe_allow_html=True)

with p4:
    st.markdown("""
<div class="card-center">
<div style="font-size:32px;">04</div>
<h4>Quantification</h4>
<p>featureCounts</p>
</div>
""", unsafe_allow_html=True)

with p5:
    st.markdown("""
<div class="card-center">
<div style="font-size:32px;">05</div>
<h4>Analysis</h4>
<p>DESeq2 + GSEA</p>
</div>
""", unsafe_allow_html=True)




# ---------- Engineering ----------

# ---------- Engineering ----------
st.write("")
st.write("")

st.markdown("## Engineering Principles")

st.markdown("<div style='margin-bottom:18px;'></div>", unsafe_allow_html=True)

e1, e2, e3 = st.columns(3, gap="large")

with e1:
    st.markdown("""
    <div class="card-center">
        <div style="font-size:48px;margin-top:6px;">🧬</div>
        <h4>Scientific</h4>
        <p>
        Reproducible statistical analysis using FastQC, DESeq2,
        GO, KEGG and GSEA.
        </p>
    </div>
    """, unsafe_allow_html=True)

with e2:
    st.markdown("""
    <div class="card-center">
        <div style="font-size:48px;margin-top:6px;">⚙️</div>
        <h4>Engineering</h4>
        <p>
        Modular Nextflow DSL2 architecture with reusable workflow
        modules and execution profiles.
        </p>
    </div>
    """, unsafe_allow_html=True)

with e3:
    st.markdown("""
    <div class="card-center">
        <div style="font-size:48px;margin-top:6px;">🗄️</div>
        <h4>Data Platform</h4>
        <p>
        SQLite-backed analytical tables supporting SQL queries
        and interactive dashboard visualization.
        </p>
    </div>
    """, unsafe_allow_html=True)

st.write("")

st.markdown("""
<div style="
background:#EAF3FF;
padding:18px 20px;
border-radius:14px;
border-left:5px solid #2563EB;
color:#1E3A8A;
font-size:16px;
font-weight:500;">
<b>Explore the complete analysis</b><br>
Use the left sidebar to navigate Quality Control, Differential Expression,
Functional Enrichment, GSEA and Downloads.
</div>
""", unsafe_allow_html=True)
