import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(page_title="NGS QC Dashboard", layout="wide")

st.title("🧬 NGS QC Dashboard")

# ========================
# 📂 Upload Section
# ========================
uploaded_file = st.file_uploader(
    "Upload QC file (CSV or TXT)",
    type=["csv", "txt"]
)

if uploaded_file is not None:

    # Detect delimiter automatically
    try:
        df = pd.read_csv(uploaded_file)
    except:
        df = pd.read_csv(uploaded_file, sep="\t")

    st.success("File uploaded successfully!")

    # ========================
    # 🧹 Clean column names
    # ========================
    df.columns = df.columns.str.strip()

    # ========================
    # 🔹 Project selection
    # ========================
    if "Project" not in df.columns:
        df["Project"] = "Default_Project"

    projects = df["Project"].unique()
    selected_project = st.selectbox("Select Project", projects)

    df_proj = df[df["Project"] == selected_project]

    st.subheader(f"📊 Project: {selected_project}")

    # ========================
    # 📈 Summary Metrics
    # ========================
    col1, col2, col3, col4 = st.columns(4)

    col1.metric("Avg Q30", f"{df_proj['% > Q30'].mean():.2f}%")
    col2.metric("Avg Duplication", f"{df_proj['% Duplication'].mean():.2f}%")
    col3.metric("Avg Adapter", f"{df_proj['% Adapter'].mean():.2f}%")
    col4.metric("Avg Mapping", f"{df_proj['Mapped_Reads(%)'].mean():.2f}%")

    # ========================
    # 🚦 QC Status Function
    # ========================
    def qc_status(row):
        if row["% Adapter"] > 50:
            return "❌ Fail"
        elif row["% Adapter"] > 20:
            return "⚠️ Warning"
        else:
            return "✅ Pass"

    df_proj["QC_Status"] = df_proj.apply(qc_status, axis=1)

    # ========================
    # 📋 Table View
    # ========================
    st.subheader("📋 Sample QC Table")
    st.dataframe(df_proj)

    # ========================
    # 📊 Visualization
    # ========================
    st.subheader("📊 QC Visualization")

    col1, col2 = st.columns(2)

    with col1:
        st.bar_chart(df_proj.set_index("Sample")["% Adapter"])

    with col2:
        st.bar_chart(df_proj.set_index("Sample")["Filtering_Rate"].str.replace("%","").astype(float))

    # ========================
    # 🧠 Insight Section
    # ========================
    st.subheader("🧠 QC Interpretation")

    avg_adapter = df_proj["% Adapter"].mean()
    avg_loss = df_proj["Filtering_Rate"].str.replace("%","").astype(float).mean()

    if avg_adapter > 50:
        st.error("High adapter contamination detected → Library prep issue")
    elif avg_adapter > 20:
        st.warning("Moderate adapter contamination → Needs optimization")
    else:
        st.success("Adapter level is acceptable")

    if avg_loss > 40:
        st.error("High read loss after trimming (>40%)")
    else:
        st.success("Read retention is good")

    # ========================
    # 🛠 Recommendation
    # ========================
    st.subheader("🛠 Recommendation")

    if avg_adapter > 50:
        st.markdown("""
        - Reduce adapter concentration
        - Improve size selection (AMPure)
        - Check adapter dimer (Bioanalyzer)
        """)

    if df_proj["% Duplication"].mean() > 25:
        st.markdown("""
        - Reduce PCR cycles
        - Increase input DNA
        """)

else:
    st.info("Please upload a QC file to begin.")
