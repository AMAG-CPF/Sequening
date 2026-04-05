# =========================================================
# IMPORTS
# =========================================================
import io
import re
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import psycopg2
import streamlit as st
from psycopg2.extras import RealDictCursor


# =========================================================
# PAGE CONFIG
# =========================================================
st.set_page_config(page_title="NGS QC App", layout="wide")
st.title("NGS QC Dashboard")


# =========================
# SECURITY
# =========================
def check_password():
    def password_entered():
        expected = st.secrets["app"]["password"]
        entered = st.session_state.get("password", "")
        st.session_state["password_correct"] = hmac.compare_digest(entered, expected)

    if st.session_state.get("password_correct", False):
        return

    st.text_input("App Password", type="password", key="password", on_change=password_entered)

    if "password" in st.session_state and st.session_state["password"] != "":
        st.error("Incorrect password")

    st.stop()


check_password()


# =========================
# DATABASE
# =========================
def get_connection():
    return psycopg2.connect(
        host=st.secrets["postgres"]["host"],
        port=st.secrets["postgres"]["port"],
        dbname=st.secrets["postgres"]["dbname"],
        user=st.secrets["postgres"]["user"],
        password=st.secrets["postgres"]["password"],
        sslmode="require",
    )


# =========================================================
# HELPER FUNCTIONS
# =========================================================
def to_numeric_safe(x):
    return pd.to_numeric(x, errors="coerce")


def clean_sample_name(x):
    return str(x).strip()


def natural_sort_key(s):
    s = str(s)
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", s)]


def smart_sort(df: pd.DataFrame) -> pd.DataFrame:
    if "Sample" in df.columns:
        return df.sort_values("Sample", key=lambda s: s.map(natural_sort_key)).reset_index(drop=True)
    return df


def format_reads(x):
    x = pd.to_numeric(x, errors="coerce")
    if pd.isna(x):
        return "NA"
    return f"{int(x):,}"


def qc_status(mapped_reads):
    x = pd.to_numeric(mapped_reads, errors="coerce")
    if pd.isna(x):
        return "Unknown"
    return "Pass" if x >= 100000 else "Maybe Low Coverage"


# =========================================================
# FILE READERS / PARSERS
# =========================================================
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


def parse_total_report_csv(uploaded_file):
    # ใช้ logic เดิมของคุณตรงนี้
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

    if "Items" in report_df.columns and "Value" in report_df.columns:
        tmp = report_df[["Items", "Value"]].copy()
    elif report_df.shape[1] >= 2:
        tmp = report_df.iloc[:, :2].copy()
        tmp.columns = ["Items", "Value"]
    else:
        return result

    for _, row in tmp.iterrows():
        key = str(row["Items"]).strip()
        value = str(row["Value"]).strip()
        if key in result:
            result[key] = value

    return result


# =========================================================
# BUILD RECORDS FOR DB INSERT
# =========================================================
def build_sample_records(df: pd.DataFrame):
    records = []
    for _, row in df.iterrows():
        records.append({
            "sample_name": row.get("Sample"),
            "raw_reads_before": pd.to_numeric(row.get("Raw_Reads_Before"), errors="coerce"),
            "duplication_percent": pd.to_numeric(row.get("% Duplication"), errors="coerce"),
            "q30_percent": pd.to_numeric(row.get("% > Q30"), errors="coerce"),
            "mb_q30_bases": pd.to_numeric(row.get("Mb Q30 bases"), errors="coerce"),
            "gc_content": pd.to_numeric(row.get("GC content"), errors="coerce"),
            "adapter_percent": pd.to_numeric(row.get("% Adapter"), errors="coerce"),
            "clean_reads_after": pd.to_numeric(row.get("Clean_Reads_After"), errors="coerce"),
            "filtering_rate": pd.to_numeric(row.get("Filtering_Rate"), errors="coerce"),
            "total_alignments": pd.to_numeric(row.get("Total_Alignments"), errors="coerce"),
            "mapped_reads": pd.to_numeric(row.get("Mapped_Reads"), errors="coerce"),
            "mapped_reads_percent": pd.to_numeric(row.get("Mapped_Reads(%)"), errors="coerce"),
            "qc_status": row.get("QC_Status"),
        })
    return records


def build_coverage_records(coverage_df: pd.DataFrame):
    if coverage_df is None or coverage_df.empty or "CHR" not in coverage_df.columns:
        return []

    records = []
    sample_cols = [c for c in coverage_df.columns if c != "CHR"]

    for _, row in coverage_df.iterrows():
        target_name = row["CHR"]
        for sample in sample_cols:
            records.append({
                "sample_name": sample,
                "target_name": target_name,
                "coverage": pd.to_numeric(row[sample], errors="coerce"),
            })

    return records


# =========================================================
# PDF REPORT
# =========================================================
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


# =========================================================
# DATABASE WRITE FUNCTIONS
# =========================================================
def insert_qc_run(run_record, sample_records, coverage_records):
    conn = get_connection()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO ngs_qc_runs (
                        project_name, flow_cell, run_date_be, software_version,
                        cycle_number, read1_length, read2_length, index1_length, index2_length,
                        total_reads_m, q30_percent, split_rate_percent, density, notes
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    RETURNING id
                    """,
                    (
                        run_record.get("project_name"),
                        run_record.get("flow_cell"),
                        run_record.get("run_date_be"),
                        run_record.get("software_version"),
                        run_record.get("cycle_number"),
                        run_record.get("read1_length"),
                        run_record.get("read2_length"),
                        run_record.get("index1_length"),
                        run_record.get("index2_length"),
                        run_record.get("total_reads_m"),
                        run_record.get("q30_percent"),
                        run_record.get("split_rate_percent"),
                        run_record.get("density"),
                        run_record.get("notes"),
                    )
                )
                run_id = cur.fetchone()[0]

                for rec in sample_records:
                    cur.execute(
                        """
                        INSERT INTO ngs_qc_samples (
                            run_id, sample_name, raw_reads_before, duplication_percent,
                            q30_percent, mb_q30_bases, gc_content, adapter_percent,
                            clean_reads_after, filtering_rate, total_alignments,
                            mapped_reads, mapped_reads_percent, qc_status
                        )
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        """,
                        (
                            run_id,
                            rec.get("sample_name"),
                            rec.get("raw_reads_before"),
                            rec.get("duplication_percent"),
                            rec.get("q30_percent"),
                            rec.get("mb_q30_bases"),
                            rec.get("gc_content"),
                            rec.get("adapter_percent"),
                            rec.get("clean_reads_after"),
                            rec.get("filtering_rate"),
                            rec.get("total_alignments"),
                            rec.get("mapped_reads"),
                            rec.get("mapped_reads_percent"),
                            rec.get("qc_status"),
                        )
                    )

                for rec in coverage_records:
                    cur.execute(
                        """
                        INSERT INTO ngs_qc_target_coverage (
                            run_id, sample_name, target_name, coverage
                        )
                        VALUES (%s, %s, %s, %s)
                        """,
                        (
                            run_id,
                            rec.get("sample_name"),
                            rec.get("target_name"),
                            rec.get("coverage"),
                        )
                    )

        return run_id
    finally:
        conn.close()


def soft_delete_qc_run(run_id: int):
    conn = get_connection()
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    UPDATE ngs_qc_runs
                    SET is_deleted = TRUE
                    WHERE id = %s
                    """,
                    (run_id,)
                )
    finally:
        conn.close()


# =========================================================
# DATABASE READ FUNCTIONS
# =========================================================
def load_qc_run_history(limit=200):
    conn = get_connection()
    try:
        query = """
            SELECT
                id,
                project_name,
                flow_cell,
                run_date_be,
                total_reads_m,
                q30_percent,
                split_rate_percent,
                density,
                created_at
            FROM ngs_qc_runs
            WHERE COALESCE(is_deleted, FALSE) = FALSE
            ORDER BY created_at DESC
            LIMIT %s
        """
        return pd.read_sql(query, conn, params=(limit,))
    finally:
        conn.close()


def load_run_info(conn, run_id):
    query = """
        SELECT
            project_name,
            flow_cell,
            run_date_be,
            software_version,
            cycle_number,
            read1_length,
            read2_length,
            index1_length,
            index2_length,
            total_reads_m,
            q30_percent,
            split_rate_percent,
            density,
            notes
        FROM ngs_qc_runs
        WHERE id = %s
    """
    df = pd.read_sql(query, conn, params=(run_id,))
    if df.empty:
        return {}

    row = df.iloc[0]
    return {
        "Project": row.get("project_name", "NA"),
        "Flow Cell": row.get("flow_cell", "NA"),
        "Run Date": row.get("run_date_be", "NA"),
        "Software Version": row.get("software_version", "NA"),
        "CycleNumber": row.get("cycle_number", "NA"),
        "Read1 Length": row.get("read1_length", "NA"),
        "Read2 Length": row.get("read2_length", "NA"),
        "Index1 Length": row.get("index1_length", "NA"),
        "Index2 Length": row.get("index2_length", "NA"),
        "TotalReads(M)": row.get("total_reads_m", "NA"),
        "Q30(%)": row.get("q30_percent", "NA"),
        "SplitRate(%)": row.get("split_rate_percent", "NA"),
        "Density(um²)": row.get("density", "NA"),
        "Notes": row.get("notes", ""),
    }


def load_sample_qc(conn, run_id):
    query = """
        SELECT
            sample_name AS sample,
            raw_reads_before AS raw_reads_before,
            duplication_percent AS duplication_percent,
            q30_percent AS q30_percent,
            mb_q30_bases AS mb_q30_bases,
            gc_content AS gc_content,
            adapter_percent AS adapter_percent,
            clean_reads_after AS clean_reads_after,
            filtering_rate AS filtering_rate,
            total_alignments AS total_alignments,
            mapped_reads AS mapped_reads,
            mapped_reads_percent AS mapped_reads_percent,
            qc_status AS qc_status
        FROM ngs_qc_samples
        WHERE run_id = %s
        ORDER BY sample_name
    """
    df = pd.read_sql(query, conn, params=(run_id,))
    if df.empty:
        return df

    return df.rename(columns={
        "sample": "Sample",
        "raw_reads_before": "Raw_Reads_Before",
        "duplication_percent": "% Duplication",
        "q30_percent": "% > Q30",
        "mb_q30_bases": "Mb Q30 bases",
        "gc_content": "GC content",
        "adapter_percent": "% Adapter",
        "clean_reads_after": "Clean_Reads_After",
        "filtering_rate": "Filtering_Rate",
        "total_alignments": "Total_Alignments",
        "mapped_reads": "Mapped_Reads",
        "mapped_reads_percent": "Mapped_Reads(%)",
        "qc_status": "QC_Status",
    })


def load_target_coverage(conn, run_id):
    query = """
        SELECT
            sample_name,
            target_name,
            coverage
        FROM ngs_qc_target_coverage
        WHERE run_id = %s
        ORDER BY sample_name, target_name
    """
    long_df = pd.read_sql(query, conn, params=(run_id,))
    if long_df.empty:
        return pd.DataFrame()

    wide_df = long_df.pivot(
        index="target_name",
        columns="sample_name",
        values="coverage"
    ).reset_index()

    wide_df = wide_df.rename(columns={"target_name": "CHR"})
    wide_df.columns.name = None
    return wide_df


# =========================================================
# RENDER FUNCTIONS
# =========================================================
def render_run_info(run_info: dict):
    if not run_info:
        st.info("No run information found.")
        return

    st.subheader("Run Information")

    r1, r2, r3 = st.columns(3)
    r1.metric("Project", run_info.get("Project", "NA"))
    r2.metric("Flow Cell", run_info.get("Flow Cell", "NA"))
    r3.metric("Run Date", run_info.get("Run Date", "NA"))

    r4, r5, r6, r7 = st.columns(4)
    r4.metric("TotalReads(M)", run_info.get("TotalReads(M)", "NA"))
    r5.metric("Q30(%)", run_info.get("Q30(%)", "NA"))
    r6.metric("SplitRate(%)", run_info.get("SplitRate(%)", "NA"))
    r7.metric("Density(um²)", run_info.get("Density(um²)", "NA"))

    with st.expander("Show full run information"):
        run_info_df = pd.DataFrame({
            "Field": list(run_info.keys()),
            "Value": list(run_info.values())
        })
        st.table(run_info_df)


def render_sample_qc(df: pd.DataFrame):
    if df is None or df.empty:
        st.warning("No sample QC data found.")
        return

    df = df.copy()
    if "QC_Status" not in df.columns and "Mapped_Reads" in df.columns:
        df["QC_Status"] = df["Mapped_Reads"].apply(qc_status)

    df = smart_sort(df)

    total_samples = len(df)
    pass_samples = int((df["QC_Status"] == "Pass").sum()) if "QC_Status" in df.columns else 0
    low_samples_n = int((df["QC_Status"] == "Maybe Low Coverage").sum()) if "QC_Status" in df.columns else 0

    st.subheader("Sample QC Summary")
    c1, c2, c3 = st.columns(3)
    c1.metric("Total Samples", total_samples)
    c2.metric("Pass", pass_samples)
    c3.metric("Maybe Low Coverage", low_samples_n)

    st.dataframe(df, use_container_width=True)


def render_coverage_module(coverage_df: pd.DataFrame):
    if coverage_df is None or coverage_df.empty:
        return

    st.subheader("Coverage Analysis")

    coverage_df = coverage_df.copy()
    coverage_df.columns = [str(c).strip() for c in coverage_df.columns]

    if "CHR" not in coverage_df.columns:
        st.warning("Coverage file/table must contain 'CHR' column.")
        return

    samples = [c for c in coverage_df.columns if c != "CHR"]
    if not samples:
        st.warning("No sample columns found in coverage data.")
        return

    c1, c2 = st.columns(2)
    with c1:
        selected_sample = st.selectbox("Select Sample", samples, key="coverage_sample")
    with c2:
        expected_targets = st.number_input(
            "Expected number of targets",
            min_value=1,
            value=len(coverage_df),
            step=1,
            key="expected_targets"
        )

    threshold = st.number_input(
        "Coverage threshold",
        min_value=0,
        value=100,
        step=1,
        key="coverage_threshold"
    )

    plot_df = coverage_df[["CHR", selected_sample]].copy()
    plot_df[selected_sample] = pd.to_numeric(plot_df[selected_sample], errors="coerce").fillna(0)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(range(len(plot_df)), plot_df[selected_sample])
    ax.axhline(threshold, linestyle="--", label=f"Threshold = {threshold}")
    ax.set_title(f"Coverage per Target ({selected_sample})")
    ax.set_xlabel("Targets")
    ax.set_ylabel("Reads")
    ax.set_xticks([])
    ax.legend()
    st.pyplot(fig)

    low_targets = plot_df[plot_df[selected_sample] < threshold].copy()

    st.subheader("Coverage Summary")
    s1, s2, s3 = st.columns(3)
    s1.metric("Expected Targets", int(expected_targets))
    s2.metric("Detected Targets", int(len(plot_df) - len(low_targets)))
    s3.metric("Low Coverage Targets", int(len(low_targets)))

    if len(low_targets) > 0:
        st.subheader("Low Coverage Targets")
        st.dataframe(low_targets, use_container_width=True)

    st.subheader("All Samples Summary")
    summary = []
    for sample in samples:
        tmp = coverage_df[["CHR", sample]].copy()
        tmp[sample] = pd.to_numeric(tmp[sample], errors="coerce").fillna(0)
        low = int((tmp[sample] < threshold).sum())
        summary.append({
            "Sample": sample,
            "Low_Targets": low,
            "Detected_Targets": int(len(tmp) - low),
            "Expected_Targets": int(expected_targets),
            "Status": "OK" if low == 0 else "Missing Targets"
        })

    summary_df = pd.DataFrame(summary)
    summary_df = summary_df.sort_values("Sample", key=lambda s: s.map(natural_sort_key)).reset_index(drop=True)
    st.dataframe(summary_df, use_container_width=True)


# =========================================================
# TABS
# ต้องมีแค่ชุดเดียว
# =========================================================
tab1, tab2 = st.tabs(["Upload & Analyze", "Browse Neon Database"])


# =========================================================
# TAB 1: Upload & Analyze
# =========================================================
with tab1:
    st.subheader("Input Information")

    h1, h2, h3 = st.columns(3)
    with h1:
        project_name_input = st.text_input("Project Name", value="NGS_QC_Project").strip()
    with h2:
        flow_cell_input = st.text_input("Flow Cell", value="").strip()
    with h3:
        run_date_input = st.text_input("Run Date (DD/MM/YY, พ.ศ.)", value="").strip()

    if not project_name_input:
        project_name_input = "NGS_QC_Project"

    date_str = datetime.now().strftime("%Y%m%d")

    st.subheader("Upload Files")
    col1, col2, col3, col4, col5 = st.columns(5)

    with col1:
        multiqc_file = st.file_uploader("MultiQC TSV", type=["tsv", "txt"], key="multiqc_file")
    with col2:
        fastp_file = st.file_uploader("FASTP TSV", type=["tsv", "txt"], key="fastp_file")
    with col3:
        flagstat_file = st.file_uploader("Flagstat TSV", type=["tsv", "txt"], key="flagstat_file")
    with col4:
        report_csv_file = st.file_uploader("Total Information CSV", type=["csv", "txt"], key="report_csv_file")
    with col5:
        coverage_file = st.file_uploader("Coverage TSV", type=["tsv", "txt"], key="coverage_file")

    df = None
    coverage_df = None
    run_info = {}

    if multiqc_file and fastp_file and flagstat_file:
        multiqc = read_table_flexible(multiqc_file)
        fastp = read_fastp_table(fastp_file)
        flagstat = read_table_flexible(flagstat_file)

        multiqc.columns = [str(c).strip() for c in multiqc.columns]
        fastp.columns = [str(c).strip() for c in fastp.columns]
        flagstat.columns = [str(c).strip() for c in flagstat.columns]

        required_multiqc = ["Sample", "% Duplication", "% > Q30", "Mb Q30 bases", "GC content", "% Adapter"]
        required_fastp = ["Sample", "Raw_Reads_Before", "Clean_Reads_After", "Filtering_Rate"]
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

        df = fastp.merge(multiqc, on="Sample", how="inner").merge(flagstat, on="Sample", how="inner")

        if "Mapped_Reads" in df.columns and "Total_Alignments" in df.columns:
            df["Mapped_Reads(%)"] = (
                pd.to_numeric(df["Mapped_Reads"], errors="coerce")
                / pd.to_numeric(df["Total_Alignments"], errors="coerce")
            ) * 100

        df["QC_Status"] = df["Mapped_Reads"].apply(qc_status)
        df = smart_sort(df)

        if report_csv_file:
            run_info = parse_total_report_csv(report_csv_file)

        if coverage_file:
            coverage_df = read_table_flexible(coverage_file)
            coverage_df.columns = [str(c).strip() for c in coverage_df.columns]

        render_run_info(run_info)
        render_sample_qc(df)

        if coverage_df is not None and not coverage_df.empty:
            render_coverage_module(coverage_df)

        st.subheader("Export QC Report")

        csv_filename = f"{project_name_input}_{date_str}_qc_report.csv"
        pdf_filename = f"{project_name_input}_{date_str}_qc_report.pdf"

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
            if col in pdf_df.columns:
                pdf_df[col] = pdf_df[col].apply(format_reads)

        for col in ["% Duplication", "% > Q30", "GC content", "% Adapter", "Mapped_Reads(%)"]:
            if col in pdf_df.columns:
                pdf_df[col] = pdf_df[col].map(lambda x: f"{to_numeric_safe(x):.2f}")

        if "Mb Q30 bases" in pdf_df.columns:
            pdf_df["Mb Q30 bases"] = pdf_df["Mb Q30 bases"].map(lambda x: f"{to_numeric_safe(x):.3f}")

        if "Filtering_Rate" in pdf_df.columns:
            pdf_df["Filtering_Rate"] = pdf_df["Filtering_Rate"].map(lambda x: f"{to_numeric_safe(x):.2f}%")

        pdf_bytes = make_pdf_report(
            pdf_df[["Sample", "Clean_Reads_After", "Mapped_Reads", "Mapped_Reads(%)", "QC_Status"]],
            project_name=project_name_input
        )

        st.download_button(
            label="Download QC Report (PDF)",
            data=pdf_bytes,
            file_name=pdf_filename,
            mime="application/pdf",
        )

        if st.button("Save Current QC Run to Neon", type="primary"):
            try:
                run_record = {
                    "project_name": project_name_input,
                    "flow_cell": flow_cell_input,
                    "run_date_be": run_date_input,
                    "software_version": run_info.get("Software Version"),
                    "cycle_number": pd.to_numeric(run_info.get("CycleNumber"), errors="coerce"),
                    "read1_length": pd.to_numeric(run_info.get("Read1 Length"), errors="coerce"),
                    "read2_length": pd.to_numeric(run_info.get("Read2 Length"), errors="coerce"),
                    "index1_length": pd.to_numeric(run_info.get("Index1 Length"), errors="coerce"),
                    "index2_length": pd.to_numeric(run_info.get("Index2 Length"), errors="coerce"),
                    "total_reads_m": pd.to_numeric(run_info.get("TotalReads(M)"), errors="coerce"),
                    "q30_percent": pd.to_numeric(run_info.get("Q30(%)"), errors="coerce"),
                    "split_rate_percent": pd.to_numeric(run_info.get("SplitRate(%)"), errors="coerce"),
                    "density": pd.to_numeric(run_info.get("Density(um²)"), errors="coerce"),
                    "notes": "",
                }

                sample_records = build_sample_records(df)
                coverage_records = build_coverage_records(coverage_df) if coverage_df is not None else []

                run_id = insert_qc_run(run_record, sample_records, coverage_records)
                st.success(f"QC run saved successfully. Run ID = {run_id}")
            except Exception as e:
                st.error(f"Failed to save QC run: {e}")


# =========================================================
# TAB 2: Browse Neon Database
# =========================================================
with tab2:
    st.subheader("Project Selection")

    try:
        run_history_df = load_qc_run_history()
    except Exception as e:
        st.error(f"Cannot connect to Neon database: {e}")
        st.stop()

    if run_history_df.empty:
        st.warning("No runs found.")
        st.stop()

    if "confirm_delete_run_id" not in st.session_state:
        st.session_state.confirm_delete_run_id = None

    project_search = st.text_input("Search Project", key="project_search_tab2").strip()

    filtered_history_df = run_history_df.copy()
    if project_search:
        filtered_history_df = filtered_history_df[
            filtered_history_df["project_name"].astype(str).str.contains(project_search, case=False, na=False)
        ]

    project_options = filtered_history_df["project_name"].dropna().astype(str).unique().tolist()

    if not project_options:
        st.warning("No matching projects found.")
        st.stop()

    selected_project = st.selectbox(
        "Select Project",
        project_options,
        key="selected_project_tab2"
    )

    project_runs = filtered_history_df[
        filtered_history_df["project_name"].astype(str) == str(selected_project)
    ].copy()

    if project_runs.empty:
        st.warning("No runs found for selected project.")
        st.stop()

    project_runs["run_label"] = project_runs.apply(
        lambda x: f"{str(x['flow_cell']) if pd.notna(x['flow_cell']) else 'NA'} | "
                  f"{str(x['run_date_be']) if pd.notna(x['run_date_be']) else 'NA'} | "
                  f"Run ID: {int(x['id'])}",
        axis=1
    )

    selected_run_label = st.selectbox(
        "Select Run",
        project_runs["run_label"].tolist(),
        key="selected_run_label_tab2"
    )

    selected_run_id = int(
        project_runs.loc[
            project_runs["run_label"] == selected_run_label, "id"
        ].iloc[0]
    )

    try:
        conn = get_connection()
        try:
            run_info = load_run_info(conn, selected_run_id)
            sample_qc_df = load_sample_qc(conn, selected_run_id)
            coverage_df_db = load_target_coverage(conn, selected_run_id)
        finally:
            conn.close()
    except Exception as e:
        st.error(f"Failed to load run data: {e}")
        st.stop()

    st.markdown("---")
    st.subheader(f"Selected Run: {selected_run_id}")

    render_run_info(run_info)

    if sample_qc_df is not None and not sample_qc_df.empty:
        render_sample_qc(sample_qc_df)
    else:
        st.info("No sample QC data found for this run.")

    if coverage_df_db is not None and not coverage_df_db.empty:
        render_coverage_module(coverage_df_db)
    else:
        st.info("No target coverage data found for this run.")

    st.markdown("---")
    st.subheader("Manage Database")

    with st.expander("🗑️ Manage Historical Runs"):
        delete_options = project_runs[["id", "project_name", "run_date_be", "flow_cell"]].copy()

        delete_options["label"] = delete_options.apply(
            lambda x: (
                f"Run ID {int(x['id'])} | "
                f"{x['project_name']} | "
                f"{x['run_date_be']} | "
                f"{x['flow_cell']}"
            ),
            axis=1
        )

        selected_delete_label = st.selectbox(
            "Select a run to delete",
            options=delete_options["label"].tolist(),
            key="delete_run_selectbox"
        )

        delete_run_id = int(
            delete_options.loc[
                delete_options["label"] == selected_delete_label, "id"
            ].iloc[0]
        )

        cdel1, cdel2 = st.columns(2)

        with cdel1:
            if st.button("Prepare Delete", type="secondary", key="prepare_delete_btn"):
                st.session_state.confirm_delete_run_id = delete_run_id

        with cdel2:
            if (
                st.session_state.confirm_delete_run_id == delete_run_id
                and st.button("Confirm Soft Delete", type="primary", key="confirm_delete_btn")
            ):
                try:
                    soft_delete_qc_run(delete_run_id)
                    st.success(f"Run ID {delete_run_id} has been deleted.")
                    st.session_state.confirm_delete_run_id = None
                    st.rerun()
                except Exception as e:
                    st.error(f"Failed to delete run: {e}")
