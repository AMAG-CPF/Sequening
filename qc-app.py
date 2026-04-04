import io
import re
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from matplotlib.ticker import FuncFormatter


st.set_page_config(page_title="NGS QC Dashboard", layout="wide")


# ==============================
# Header
# ==============================
st.title("🧬 NGS QC Dashboard")

col_a, col_b, col_c = st.columns(3)

with col_a:
    project_name_input = st.text_input(
        "Project Name (default if not in file)",
        value="NGS_QC_Project"
    ).strip()

with col_b:
    flow_cell_input = st.text_input(
        "Flow Cell",
        value=""
    ).strip()

with col_c:
    run_date_input = st.text_input(
        "Run Date (DD/MM/YY, พ.ศ.)",
        value=""
    ).strip()

if not project_name_input:
    project_name_input = "NGS_QC_Project"

date_str = datetime.now().strftime("%Y%m%d")

# ==============================
# Upload Section
# ==============================
col1, col2, col3, col4 = st.columns(4)

with col1:
    multiqc_file = st.file_uploader("Upload MultiQC TSV", type=["tsv", "txt"])

with col2:
    fastp_file = st.file_uploader("Upload FASTP TSV", type=["tsv", "txt"])

with col3:
    flagstat_file = st.file_uploader("Upload Flagstat TSV", type=["tsv", "txt"])

with col4:
    report_csv_file = st.file_uploader("Upload Total Information CSV", type=["csv", "txt"])

# ==============================
# Helper Functions
# ==============================
import ast

def parse_total_report_csv(uploaded_file):
    uploaded_file.seek(0)

    result = {
        "Software Version": "NA",
        "CycleNumber": "NA",
        "Read1 Length": "NA",
        "Read2 Length": "NA",
        "Index1 Length": "NA",
        "Index2 Length": "NA",
        "Read1 Dark Length": "NA",
        "Read2 Dark Length": "NA",
        "TotalReads(M)": "NA",
        "Q30(%)": "NA",
        "SplitRate(%)": "NA",
        "Density(um²)": "NA",
        "R1 Phasing": "NA",
        "R1 Prephasing": "NA",
        "R2 Phasing": "NA",
        "R2 Prephasing": "NA",
        "MaxOffsetX": "NA",
        "MaxOffsetY": "NA",
    }

    try:
        report_df = pd.read_csv(uploaded_file, dtype=str)
    except Exception:
        uploaded_file.seek(0)
        report_df = pd.read_csv(uploaded_file, sep="\t", dtype=str)

    report_df.columns = [str(c).strip() for c in report_df.columns]

    # รองรับทั้งแบบมี header "Items,Value"
    if "Items" in report_df.columns and "Value" in report_df.columns:
        for _, row in report_df.iterrows():
            key = str(row["Items"]).strip()
            value = str(row["Value"]).strip()
            if key in result:
                result[key] = value

    # รองรับแบบไม่มี header หรือมีแค่ 2 คอลัมน์
    elif report_df.shape[1] >= 2:
        tmp = report_df.iloc[:, :2].copy()
        tmp.columns = ["Items", "Value"]
        for _, row in tmp.iterrows():
            key = str(row["Items"]).strip()
            value = str(row["Value"]).strip()
            if key in result:
                result[key] = value

    return result

def clean_sample_name(sample: str) -> str:
    sample = str(sample).strip()
    sample = re.sub(r"\.trimmed$", "", sample)
    sample = re.sub(r"\.flagstats$", "", sample)
    return sample


def to_numeric_safe(x):
    if pd.isna(x):
        return 0.0
    x = str(x).replace(",", "").replace("%", "").strip()
    try:
        return float(x)
    except Exception:
        return 0.0


def format_reads(x) -> str:
    x = to_numeric_safe(x)
    if x >= 1_000_000:
        return f"{x / 1_000_000:.1f}M"
    if x >= 1_000:
        return f"{x / 1_000:.1f}k"
    return str(int(round(x)))


def qc_status(mapped_reads) -> str:
    mapped_reads = to_numeric_safe(mapped_reads)
    return "Maybe Low Coverage" if mapped_reads < 10000 else "Pass"


def natural_sort_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def smart_sort(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(
        by="Sample",
        key=lambda s: s.map(natural_sort_key)
    ).reset_index(drop=True)


def read_table_flexible(uploaded_file) -> pd.DataFrame:
    uploaded_file.seek(0)
    try:
        df = pd.read_csv(uploaded_file, sep="\t", dtype=str)
        if df.shape[1] == 1:
            uploaded_file.seek(0)
            df = pd.read_csv(uploaded_file, sep=r"\s+", engine="python", dtype=str)
    except Exception:
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file, sep=r"\s+", engine="python", dtype=str)
    return df


def read_fastp_table(uploaded_file) -> pd.DataFrame:
    uploaded_file.seek(0)
    fastp = pd.read_csv(uploaded_file, sep="\t", header=None, dtype=str)

    if fastp.shape[1] == 1:
        uploaded_file.seek(0)
        fastp = pd.read_csv(uploaded_file, sep=r"\s+", header=None, engine="python", dtype=str)

    if fastp.shape[1] == 4:
        fastp.columns = ["Sample", "Raw_Reads_Before", "Clean_Reads_After", "Filtering_Rate"]
    elif fastp.shape[1] == 3:
        fastp.columns = ["Sample", "Raw_Reads_Before", "Clean_Reads_After"]
        fastp["Filtering_Rate"] = 0
    else:
        st.error(f"FASTP file has {fastp.shape[1]} columns, expected 3 or 4 columns.")
        st.stop()

    return fastp[["Sample", "Raw_Reads_Before", "Clean_Reads_After", "Filtering_Rate"]]


def make_pdf_report(df_report: pd.DataFrame, project_name: str) -> io.BytesIO:
    buffer = io.BytesIO()

    n_rows = len(df_report)
    fig_height = max(3, 1 + 0.45 * (n_rows + 2))

    fig, ax = plt.subplots(figsize=(12, fig_height))
    ax.axis("off")
    ax.set_title(f"{project_name} - QC Summary Report", fontsize=14, fontweight="bold", pad=16)

    table = ax.table(
        cellText=df_report.values,
        colLabels=df_report.columns,
        loc="center",
        cellLoc="center",
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)

    status_col = list(df_report.columns).index("QC_Status") if "QC_Status" in df_report.columns else None

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold")
        elif status_col is not None and col == status_col:
            if str(df_report.iloc[row - 1, status_col]) == "Maybe Low Coverage":
                cell.set_facecolor("#f8d7da")
            else:
                cell.set_facecolor("#d1e7dd")

    fig.tight_layout()
    fig.savefig(buffer, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return buffer


def highlight_status(row):
    if row["QC_Status"] == "Maybe Low Coverage":
        return ["background-color: #ffe6e6"] * len(row)
    return ["background-color: #e6ffe6"] * len(row)


# ==============================
# Main Logic
# ==============================
if multiqc_file and fastp_file and flagstat_file:
    multiqc = read_table_flexible(multiqc_file)
    fastp = read_fastp_table(fastp_file)
    flagstat = read_table_flexible(flagstat_file)

    multiqc.columns = [str(c).strip() for c in multiqc.columns]
    fastp.columns = [str(c).strip() for c in fastp.columns]
    flagstat.columns = [str(c).strip() for c in flagstat.columns]

    required_multiqc = [
        "Sample", "% Duplication", "% > Q30", "Mb Q30 bases",
        "GC content", "% Adapter"
    ]
    required_fastp = [
        "Sample", "Raw_Reads_Before", "Clean_Reads_After", "Filtering_Rate"
    ]
    required_flagstat = ["Sample", "Total_Alignments", "Mapped_Reads"]

    for col in required_multiqc:
        if col not in multiqc.columns:
            st.error(f"Missing column in MultiQC file: {col}")
            st.stop()

    for col in required_fastp:
        if col not in fastp.columns:
            st.error(f"Missing column in FASTP file: {col}")
            st.stop()

    for col in required_flagstat:
        if col not in flagstat.columns:
            st.error(f"Missing column in Flagstat file: {col}")
            st.stop()

    multiqc["Sample"] = multiqc["Sample"].apply(clean_sample_name)
    fastp["Sample"] = fastp["Sample"].apply(clean_sample_name)
    flagstat["Sample"] = flagstat["Sample"].apply(clean_sample_name)

    multiqc = multiqc[required_multiqc].copy()
    fastp = fastp[required_fastp].copy()
    flagstat = flagstat[required_flagstat].copy()

    for col in ["% Duplication", "% > Q30", "Mb Q30 bases", "GC content", "% Adapter"]:
        multiqc[col] = multiqc[col].apply(to_numeric_safe)

    for col in ["Raw_Reads_Before", "Clean_Reads_After", "Filtering_Rate"]:
        fastp[col] = fastp[col].apply(to_numeric_safe)

    for col in ["Total_Alignments", "Mapped_Reads"]:
        flagstat[col] = flagstat[col].apply(to_numeric_safe)

    df = multiqc.merge(fastp, on="Sample", how="inner")
    df = df.merge(flagstat, on="Sample", how="inner")

    if df.empty:
        st.error("No matching samples found after merging. Please check sample names in the 3 files.")
        st.stop()

    df["Mapped_Reads(%)"] = df.apply(
        lambda row: (row["Mapped_Reads"] / row["Total_Alignments"] * 100)
        if row["Total_Alignments"] > 0 else 0,
        axis=1
    )

    # ==============================
    # Project handling
    # ==============================
    if "Project" not in df.columns:
        df["Project"] = project_name_input

    projects = sorted(df["Project"].dropna().unique(), key=natural_sort_key)

    selected_project = st.selectbox("Select Project", projects)
    df = df[df["Project"] == selected_project].copy()

    st.markdown(f"### 📁 Project: {selected_project}")

# ==============================
# Run Information
# ==============================
st.subheader("Run Information")

run_info = {
    "Project": selected_project,
    "Flow Cell": flow_cell_input if flow_cell_input else "NA",
    "Run Date": run_date_input if run_date_input else "NA",
    "Software Version": "NA",
    "CycleNumber": "NA",
    "Read1 Length": "NA",
    "Read2 Length": "NA",
    "Index1 Length": "NA",
    "Index2 Length": "NA",
    "Read1 Dark Length": "NA",
    "Read2 Dark Length": "NA",
    "TotalReads(M)": "NA",
    "Q30(%)": "NA",
    "SplitRate(%)": "NA",
    "Density(um²)": "NA",
    "R1 Phasing": "NA",
    "R1 Prephasing": "NA",
    "R2 Phasing": "NA",
    "R2 Prephasing": "NA",
    "MaxOffsetX": "NA",
    "MaxOffsetY": "NA",
}

if report_csv_file is not None:
    parsed_info = parse_total_report_csv(report_csv_file)
    run_info.update(parsed_info)

r1, r2, r3 = st.columns(3)
r1.metric("Project", run_info["Project"])
r2.metric("Flow Cell", run_info["Flow Cell"])
r3.metric("Run Date", run_info["Run Date"])

r4, r5, r6, r7 = st.columns(4)
r4.metric("TotalReads(M)", run_info["TotalReads(M)"])
r5.metric("Q30(%)", run_info["Q30(%)"])
r6.metric("SplitRate(%)", run_info["SplitRate(%)"])
r7.metric("Density(um²)", run_info["Density(um²)"])

with st.expander("Show full run information"):
    run_info_df = pd.DataFrame({
        "Field": list(run_info.keys()),
        "Value": list(run_info.values())
    })
    st.table(run_info_df)

    # ==============================
    # QC and arrangement
    # ==============================
    df["QC_Status"] = df["Mapped_Reads"].apply(qc_status)

    df = df[
        [
            "Sample",
            "Raw_Reads_Before",
            "% Duplication",
            "% > Q30",
            "Mb Q30 bases",
            "GC content",
            "% Adapter",
            "Clean_Reads_After",
            "Filtering_Rate",
            "Total_Alignments",
            "Mapped_Reads",
            "Mapped_Reads(%)",
            "QC_Status",
            "Project",
        ]
    ]

    df = smart_sort(df)

    # ==============================
    # Summary
    # ==============================
    total_samples = len(df)
    pass_samples = (df["QC_Status"] == "Pass").sum()
    low_samples_n = (df["QC_Status"] == "Maybe Low Coverage").sum()

    st.subheader("Summary")
    c1, c2, c3 = st.columns(3)
    c1.metric("Total Samples", int(total_samples))
    c2.metric("Pass", int(pass_samples))
    c3.metric("Maybe Low Coverage", int(low_samples_n))

    # ==============================
    # QC Highlight Table
    # ==============================
    st.subheader("QC Status Table")

    df_display = df.copy()

    percent_cols = ["% Duplication", "% > Q30", "GC content", "% Adapter", "Mapped_Reads(%)"]
    for col in percent_cols:
        df_display[col] = df_display[col].map(lambda x: f"{x:.2f}")

    df_display["Mb Q30 bases"] = df_display["Mb Q30 bases"].map(lambda x: f"{x:.3f}")
    df_display["Filtering_Rate"] = df_display["Filtering_Rate"].map(lambda x: f"{x:.2f}%")

    read_cols = ["Raw_Reads_Before", "Clean_Reads_After", "Total_Alignments", "Mapped_Reads"]
    for col in read_cols:
        df_display[col] = df_display[col].apply(format_reads)

    st.dataframe(
        df_display.style.apply(highlight_status, axis=1),
        use_container_width=True
    )

    # ==============================
    # Plot with threshold line
    # ==============================
    st.subheader("Reads Distribution")

    fig, ax = plt.subplots(figsize=(12, 5))

    plot_df = df.copy()
    x = range(len(plot_df))

    ax.bar(
        [i - 0.2 for i in x],
        plot_df["Clean_Reads_After"],
        width=0.4,
        label="Clean Reads"
    )
    ax.bar(
        [i + 0.2 for i in x],
        plot_df["Mapped_Reads"],
        width=0.4,
        label="Mapped Reads"
    )

    ax.axhline(10000, linestyle="--", linewidth=1.5, label="Threshold = 10,000")

    ax.set_xticks(list(x))
    ax.set_xticklabels(plot_df["Sample"], rotation=45, ha="right")
    ax.set_ylabel("Reads")
    ax.set_title("Clean Reads and Mapped Reads by Sample")
    ax.legend()

    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda y, _: f"{y/1000:.0f}k" if y >= 1000 else f"{int(y)}")
    )

    fig.tight_layout()
    st.pyplot(fig)

    # ==============================
    # Flagged samples
    # ==============================
    low_samples = df[df["QC_Status"] == "Maybe Low Coverage"].copy()

    if len(low_samples) > 0:
        st.subheader("Flagged Samples")
        low_display = low_samples[["Sample", "Clean_Reads_After", "Mapped_Reads", "QC_Status"]].copy()
        low_display["Clean_Reads_After"] = low_display["Clean_Reads_After"].apply(format_reads)
        low_display["Mapped_Reads"] = low_display["Mapped_Reads"].apply(format_reads)
        st.dataframe(low_display, use_container_width=True)
    else:
        st.success("All samples passed the minimum mapped reads threshold.")

    # ==============================
    # Export
    # ==============================
    st.subheader("Export QC Report")

    csv_filename = f"{selected_project}_{date_str}_qc_report.csv"
    pdf_filename = f"{selected_project}_{date_str}_qc_report.pdf"

    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)

    st.download_button(
        label="Download QC Report (CSV)",
        data=csv_buffer.getvalue(),
        file_name=csv_filename,
        mime="text/csv",
    )

    pdf_df = df.copy()
    for col in ["Raw_Reads_Before", "Clean_Reads_After", "Total_Alignments", "Mapped_Reads"]:
        pdf_df[col] = pdf_df[col].apply(format_reads)

    for col in ["% Duplication", "% > Q30", "GC content", "% Adapter", "Mapped_Reads(%)"]:
        pdf_df[col] = pdf_df[col].map(lambda x: f"{x:.2f}")

    pdf_df["Mb Q30 bases"] = pdf_df["Mb Q30 bases"].map(lambda x: f"{x:.3f}")
    pdf_df["Filtering_Rate"] = pdf_df["Filtering_Rate"].map(lambda x: f"{x:.2f}%")

    pdf_bytes = make_pdf_report(
        pdf_df[["Sample", "Clean_Reads_After", "Mapped_Reads", "Mapped_Reads(%)", "QC_Status"]],
        project_name=selected_project
    )

    st.download_button(
        label="Download QC Report (PDF)",
        data=pdf_bytes,
        file_name=pdf_filename,
        mime="application/pdf",
    )
