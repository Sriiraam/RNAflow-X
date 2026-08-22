import streamlit as st
import pandas as pd
import plotly.express as px

st.set_page_config(layout="wide")

# -----------------------------
# Custom CSS
# -----------------------------
st.markdown("""
<style>
.hero-qc{
    padding:34px;
    border-radius:22px;
    background:linear-gradient(135deg,#021B49 0%,#0B4CC2 100%);
    color:white;
    box-shadow:0 10px 30px rgba(0,70,180,.25);
}
.hero-qc h1{
    font-size:46px;
    margin:0;
    font-weight:800;
}
.hero-qc h3{
    color:#67E8F9;
    margin:8px 0 14px 0;
    font-weight:600;
}
.hero-qc p{
    color:#E5F3FF;
    font-size:16px;
    line-height:1.7;
}

.metric-card{
    background:white;
    border-radius:18px;
    padding:18px;
    text-align:center;
    border:1px solid #E5E7EB;
    box-shadow:0 4px 14px rgba(0,0,0,.06);
}
.metric-label{
    color:#64748B;
    font-size:14px;
}
.metric-value{
    color:#0F172A;
    font-size:30px;
    font-weight:800;
}

.section-title{
    font-size:28px;
    font-weight:800;
    margin:22px 0 12px 0;
}
</style>
""", unsafe_allow_html=True)

# -----------------------------
# Hero Banner
# -----------------------------
st.markdown("""
<div class="hero-qc">
    <h1>Quality Control Dashboard</h1>
    <h3>FastQC • MultiQC • Sequencing Quality Assessment</h3>
    <p>
    This dashboard summarizes sequencing quality before transcript
    quantification. Metrics are aggregated from MultiQC and FastQC.
    </p>
</div>
""", unsafe_allow_html=True)

st.write("")

# -----------------------------
# QC Dataset
# -----------------------------
qc = pd.DataFrame({
    "Sample":[
        "control_rep1",
        "control_rep2",
        "pfos50_rep1",
        "pfos50_rep2"
    ],
    "GC %":[48.4,43.5,51.3,50.5],
    "Read Length":["36 bp"]*4,
    "Mean Quality":[37.8,37.6,38.1,37.9]
})

# -----------------------------
# Sample Summary
# -----------------------------
st.markdown('<div class="section-title">📊 Sample Quality Summary</div>',
            unsafe_allow_html=True)

styled_qc = (
    qc.style
    .set_table_styles([
        {"selector":"th","props":[("text-align","center"),("font-size","15px")]},
        {"selector":"td","props":[("text-align","center"),("font-size","14px")]},
        {"selector":"th:nth-child(1)","props":[("color","#2563EB"),("font-weight","700")]},
        {"selector":"th:nth-child(2)","props":[("color","#16A34A"),("font-weight","700")]},
        {"selector":"th:nth-child(3)","props":[("color","#7C3AED"),("font-weight","700")]},
        {"selector":"th:nth-child(4)","props":[("color","#EA580C"),("font-weight","700")]},
    ])
    .set_properties(subset=["GC %"], **{
        "background-color":"#ECFDF3",
        "color":"#15803D",
        "font-weight":"700"
    })
    .set_properties(subset=["Read Length"], **{
        "background-color":"#F5F3FF",
        "color":"#6D28D9",
        "font-weight":"700"
    })
    .set_properties(subset=["Mean Quality"], **{
        "background-color":"#FFF7ED",
        "color":"#C2410C",
        "font-weight":"700"
    })
)

st.dataframe(styled_qc, width="stretch", hide_index=True)

st.write("")

# -----------------------------
# KPI Cards
# -----------------------------
c1,c2,c3,c4 = st.columns(4)

cards = [
    ("Samples Passed","4 / 4"),
    ("Average GC","48.4%"),
    ("Read Length","36 bp"),
    ("Mean Q Score","37.9")
]

for col,(label,value) in zip([c1,c2,c3,c4],cards):
    with col:
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">{label}</div>
            <div class="metric-value">{value}</div>
        </div>
        """, unsafe_allow_html=True)

st.write("")

# -----------------------------
# GC Plot + Interpretation
# -----------------------------
left,right = st.columns([1.1,1])

with left:
    st.markdown("### GC Content Across Samples")

    fig = px.bar(
        qc,
        x="Sample",
        y="GC %",
        color="GC %",
        color_continuous_scale="Tealgrn"
    )

    fig.update_layout(
        template="plotly_white",
        height=360,
        coloraxis_showscale=False,
        margin=dict(l=10,r=10,t=30,b=10)
    )

    st.plotly_chart(fig, width="stretch")

with right:
    st.markdown("### QC Interpretation")

    st.success("""
**Overall sequencing quality is excellent**

- All 4 libraries passed FastQC
- Mean quality exceeds Q37
- GC distribution is biologically consistent
- No adapter contamination detected
""")

    st.info("""
**Why this matters**

High-quality sequencing improves transcript quantification,
reduces false positives, and increases confidence in GO,
KEGG and GSEA analyses.
""")

st.write("")

# -----------------------------
# MultiQC Status
# -----------------------------
st.markdown('<div class="section-title">🧬 MultiQC Module Summary</div>',
            unsafe_allow_html=True)

status = pd.DataFrame({
    "Module":[
        "Per Base Quality",
        "GC Content",
        "Adapter Content",
        "N Content",
        "Sequence Quality",
        "Length Distribution"
    ],
    "Status":["✓ PASS"]*6
})

styled_status = (
    status.style
    .set_table_styles([
        {"selector":"th","props":[("text-align","center"),("font-size","15px")]},
        {"selector":"td","props":[("text-align","center"),("font-size","14px")]},
        {"selector":"th:nth-child(1)","props":[("color","#2563EB"),("font-weight","700")]},
        {"selector":"th:nth-child(2)","props":[("color","#16A34A"),("font-weight","700")]},
    ])
    .set_properties(subset=["Status"], **{
        "background-color":"#ECFDF3",
        "color":"#15803D",
        "font-weight":"800"
    })
)

st.dataframe(styled_status, width="stretch", hide_index=True)

st.caption("Source: MultiQC v1.18 • FastQC v0.12.1")
