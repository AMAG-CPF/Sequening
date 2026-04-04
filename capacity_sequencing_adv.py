import streamlit as st
import plotly.graph_objects as go
import pandas as pd

st.set_page_config(page_title="NGS Capacity Planner", layout="wide")

# -------------------------
# CONSTANTS
# -------------------------
READ_LENGTH = 300  # PE150

KIT_OPTIONS = {
    "25M (7.5 Gb)": 7.5e9,
    "60M (18 Gb)": 18e9
}

# -------------------------
# UI
# -------------------------
st.title("🧬 NGS Capacity Planning Tool (Salus Nimbo)")

col1, col2 = st.columns(2)

with col1:
    kit = st.radio("Select Reagent Kit", list(KIT_OPTIONS.keys()))
    coverage = st.slider("Coverage (X)", 5, 100, 20)

with col2:
    parentage_samples = st.number_input("Parentage samples", 0, 10000, 300)
    wssv_samples = st.number_input("WSSV samples", 0, 10000, 200)

# Panel sizes
PARENTAGE_SIZE = 30645
WSSV_SIZE = 10603

# -------------------------
# CALCULATION
# -------------------------
capacity = KIT_OPTIONS[kit]

parentage_bp = parentage_samples * PARENTAGE_SIZE * coverage
wssv_bp = wssv_samples * WSSV_SIZE * coverage

total_bp = parentage_bp + wssv_bp

# Reads
expected_reads = capacity / READ_LENGTH
required_reads = total_bp / READ_LENGTH

# Metrics
usage_percent = (total_bp / capacity) * 100
reads_usage_percent = (required_reads / expected_reads) * 100

remaining_bp = capacity - total_bp
remaining_reads = expected_reads - required_reads

total_samples = parentage_samples + wssv_samples
reads_per_sample = required_reads / total_samples if total_samples > 0 else 0

# Max samples estimation (based on current mix ratio)
if total_samples > 0:
    avg_bp_per_sample = total_bp / total_samples
    max_samples = capacity / avg_bp_per_sample
else:
    max_samples = 0

# -------------------------
# METRICS DISPLAY
# -------------------------
st.subheader("📊 Summary")

col1, col2, col3, col4 = st.columns(4)

col1.metric("Total Required (Gb)", f"{total_bp/1e9:.2f}")
col2.metric("Capacity (Gb)", f"{capacity/1e9:.2f}")
col3.metric("Usage (%)", f"{usage_percent:.1f}")
col4.metric("Remaining (Gb)", f"{remaining_bp/1e9:.2f}")

col1, col2, col3, col4 = st.columns(4)

col1.metric("Reads Required (M)", f"{required_reads/1e6:.2f}")
col2.metric("Reads Capacity (M)", f"{expected_reads/1e6:.2f}")
col3.metric("Reads / Sample", f"{reads_per_sample:.0f}")
col4.metric("Max Samples/run", f"{int(max_samples)}")

# -------------------------
# WARNING
# -------------------------
if total_bp > capacity:
    st.error("⚠️ Exceeds sequencing capacity!")
else:
    st.success("✅ Within capacity")

# -------------------------
# PLOT (Gb)
# -------------------------
fig = go.Figure()

fig.add_trace(go.Bar(
    x=["Parentage", "WSSV"],
    y=[parentage_bp/1e9, wssv_bp/1e9],
    name="Usage (Gb)"
))

fig.add_hline(y=capacity/1e9)

fig.update_layout(
    title="Sequencing Usage (Gb)",
    yaxis_title="Gb"
)

st.plotly_chart(fig, use_container_width=True)

# -------------------------
# PLOT (Reads)
# -------------------------
fig2 = go.Figure()

fig2.add_trace(go.Bar(
    x=["Required Reads"],
    y=[required_reads/1e6],
    name="Required"
))

fig2.add_trace(go.Bar(
    x=["Capacity Reads"],
    y=[expected_reads/1e6],
    name="Capacity"
))

fig2.update_layout(
    title="Reads Usage (Million)",
    yaxis_title="Reads (M)"
)

st.plotly_chart(fig2, use_container_width=True)

# -------------------------
# EXPORT
# -------------------------

st.subheader("📁 Export Report")

report_df = pd.DataFrame([{
    "Selected_Reagent_Kit": kit,
    "Coverage_X": coverage,
    "Parentage_Samples": parentage_samples,
    "WSSV_Samples": wssv_samples,
    "Total_Samples": total_samples,
    "Total_Gb_Required": total_bp/1e9,
    "Capacity_Gb": capacity/1e9,
    "Usage_Percent": usage_percent,
    "Remaining_Gb": remaining_bp/1e9,
    "Reads_Required_M": required_reads/1e6,
    "Reads_Capacity_M": expected_reads/1e6,
    "Reads_per_Sample": reads_per_sample,
    "Max_Samples": max_samples
}])

st.dataframe(report_df)

st.download_button(
    "Download CSV",
    report_df.to_csv(index=False),
    file_name="NGS_planning_report.csv",
    mime="text/csv"
)
