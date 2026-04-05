import io
import re
import numpy as np
import hmac
from datetime import datetime
import psycopg2
from psycopg2.extras import RealDictCursor
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from matplotlib.ticker import FuncFormatter
from sqlalchemy import create_engine, text


# =========================================================
# Streamlit config
# =========================================================
st.set_page_config(page_title="NGS QC Dashboard", layout="wide")
st.title("🧬 NGS QC Dashboard")


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

def insert_qc_run(run_record: dict, sample_records: list[dict], coverage_records: list[dict]) -> int:
    conn = get_connection()
    try:
        run_record = normalize_record(run_record)
        sample_records = [normalize_record(rec) for rec in sample_records]
        coverage_records = [normalize_record(rec) for rec in coverage_records]

        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO ngs_qc_runs (
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
                    )
                    VALUES (
                        %(project_name)s,
                        %(flow_cell)s,
                        %(run_date_be)s,
                        %(software_version)s,
                        %(cycle_number)s,
                        %(read1_length)s,
                        %(read2_length)s,
                        %(index1_length)s,
                        %(index2_length)s,
                        %(total_reads_m)s,
                        %(q30_percent)s,
                        %(split_rate_percent)s,
                        %(density)s,
                        %(notes)s
                    )
                    RETURNING id
                    """,
                    run_record
                )
                run_id = cur.fetchone()[0]

                for rec in sample_records:
                    rec["run_id"] = int(run_id)
                    rec = normalize_record(rec)
                    cur.execute(
                        """
                        INSERT INTO ngs_qc_samples (
                            run_id,
                            sample_name,
                            raw_reads_before,
                            duplication_percent,
                            q30_percent,
                            mb_q30_bases,
                            gc_content,
                            adapter_percent,
                            clean_reads_after,
                            filtering_rate,
                            total_alignments,
                            mapped_reads,
                            mapped_reads_percent,
                            qc_status
                        )
                        VALUES (
                            %(run_id)s,
                            %(sample_name)s,
                            %(raw_reads_before)s,
                            %(duplication_percent)s,
                            %(q30_percent)s,
                            %(mb_q30_bases)s,
                            %(gc_content)s,
                            %(adapter_percent)s,
                            %(clean_reads_after)s,
                            %(filtering_rate)s,
                            %(total_alignments)s,
                            %(mapped_reads)s,
                            %(mapped_reads_percent)s,
                            %(qc_status)s
                        )
                        """,
                        rec
                    )

                for rec in coverage_records:
                    rec["run_id"] = int(run_id)
                    rec = normalize_record(rec)
                    cur.execute(
                        """
                        INSERT INTO ngs_qc_target_coverage (
                            run_id,
                            sample_name,
                            target_name,
                            coverage
                        )
                        VALUES (
                            %(run_id)s,
                            %(sample_name)s,
                            %(target_name)s,
                            %(coverage)s
                        )
                        """,
                        rec
                    )

        return run_id
    finally:
        conn.close()

def build_sample_records(df: pd.DataFrame) -> list[dict]:
    out = []
    for _, row in df.iterrows():
        out.append({
            "sample_name": row["Sample"],
            "raw_reads_before": float(row["Raw_Reads_Before"]),
            "duplication_percent": float(row["% Duplication"]),
            "q30_percent": float(row["% > Q30"]),
            "mb_q30_bases": float(row["Mb Q30 bases"]),
            "gc_content": float(row["GC content"]),
            "adapter_percent": float(row["% Adapter"]),
            "clean_reads_after": float(row["Clean_Reads_After"]),
            "filtering_rate": float(row["Filtering_Rate"]),
            "total_alignments": float(row["Total_Alignments"]),
            "mapped_reads": float(row["Mapped_Reads"]),
            "mapped_reads_percent": float(row["Mapped_Reads(%)"]),
            "qc_status": row["QC_Status"],
        })
    return out

def build_coverage_records(coverage_df: pd.DataFrame) -> list[dict]:
    out = []
    samples = [c for c in coverage_df.columns if c != "CHR"]

    for _, row in coverage_df.iterrows():
        target_name = row["CHR"]
        for sample in samples:
            out.append({
                "sample_name": sample,
                "target_name": target_name,
                "coverage": float(pd.to_numeric(row[sample], errors="coerce") or 0),
            })
    return out

def load_qc_run_history(limit=200):
    conn = get_connection()
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(
                """
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
                """,
                (limit,)
            )
            rows = cur.fetchall()
        return pd.DataFrame(rows)
    finally:
        conn.close()
def load_qc_samples(run_id: int):
    conn = get_connection()
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(
                """
                SELECT
                    sample_name AS "Sample",
                    raw_reads_before AS "Raw_Reads_Before",
                    duplication_percent AS "% Duplication",
                    q30_percent AS "% > Q30",
                    mb_q30_bases AS "Mb Q30 bases",
                    gc_content AS "GC content",
                    adapter_percent AS "% Adapter",
                    clean_reads_after AS "Clean_Reads_After",
                    filtering_rate AS "Filtering_Rate",
                    total_alignments AS "Total_Alignments",
                    mapped_reads AS "Mapped_Reads",
                    mapped_reads_percent AS "Mapped_Reads(%)",
                    qc_status AS "QC_Status"
                FROM ngs_qc_samples
                WHERE run_id = %s
                ORDER BY sample_name
                """,
                (run_id,)
            )
            rows = cur.fetchall()
        return pd.DataFrame(rows)
    finally:
        conn.close()

def load_qc_coverage(run_id: int):
    conn = get_connection()
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(
                """
                SELECT sample_name, target_name, coverage
                FROM ngs_qc_target_coverage
                WHERE run_id = %s
                ORDER BY target_name, sample_name
                """,
                (run_id,)
            )
            rows = cur.fetchall()

        long_df = pd.DataFrame(rows)
        if long_df.empty:
            return pd.DataFrame()

        wide_df = long_df.pivot(
            index="target_name",
            columns="sample_name",
            values="coverage"
        ).reset_index()

        wide_df = wide_df.rename(columns={"target_name": "CHR"})
        return wide_df
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
# Helper functions
# =========================================================

def to_python_type(value):
    if pd.isna(value):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value

def normalize_record(record: dict) -> dict:
    return {k: to_python_type(v) for k, v in record.items()}

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


def highlight_status(row):
    if row["QC_Status"] == "Maybe Low Coverage":
        return ["background-color: #ffe6e6"] * len(row)
    return ["background-color: #e6ffe6"] * len(row)


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
# Database functions
# =========================================================
def load_projects(conn):
    query = """
        SELECT DISTINCT project_name
        FROM ngs_qc_runs
        WHERE project_name IS NOT NULL
          AND COALESCE(is_deleted, FALSE) = FALSE
        ORDER BY project_name
    """
    return pd.read_sql(query, conn)


def load_runs_by_project(conn, project_name):
    query = """
        SELECT
            id,
            flow_cell,
            run_date_be,
            software_version,
            cycle_number,
            total_reads_m,
            q30_percent,
            split_rate_percent,
            density,
            notes,
            created_at
        FROM ngs_qc_runs
        WHERE project_name = %s
          AND COALESCE(is_deleted, FALSE) = FALSE
        ORDER BY id DESC
    """
    return pd.read_sql(query, conn, params=(project_name,))

def load_run_info(conn, run_id):
    query = """
        SELECT
            id,
            created_at,
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
# Reusable render functions
# =========================================================
def render_run_info(run_info: dict):
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
    if df.empty:
        st.warning("No sample QC data found.")
        return

    df = df.copy()
    if "QC_Status" not in df.columns and "Mapped_Reads" in df.columns:
        df["QC_Status"] = df["Mapped_Reads"].apply(qc_status)

    if "Sample" in df.columns:
        df = smart_sort(df)

    total_samples = len(df)
    pass_samples = (df["QC_Status"] == "Pass").sum()
    low_samples_n = (df["QC_Status"] == "Maybe Low Coverage").sum()

    st.subheader("Summary")
    c1, c2, c3 = st.columns(3)
    c1.metric("Total Samples", int(total_samples))
    c2.metric("Pass", int(pass_samples))
    c3.metric("Maybe Low Coverage", int(low_samples_n))

    st.subheader("QC Status Table")

    df_display = df.copy()

    percent_cols = [c for c in ["% Duplication", "% > Q30", "GC content", "% Adapter", "Mapped_Reads(%)"] if c in df_display.columns]
    for col in percent_cols:
        df_display[col] = df_display[col].map(lambda x: f"{to_numeric_safe(x):.2f}")

    if "Mb Q30 bases" in df_display.columns:
        df_display["Mb Q30 bases"] = df_display["Mb Q30 bases"].map(lambda x: f"{to_numeric_safe(x):.3f}")

    if "Filtering_Rate" in df_display.columns:
        df_display["Filtering_Rate"] = df_display["Filtering_Rate"].map(lambda x: f"{to_numeric_safe(x):.2f}%")

    read_cols = [c for c in ["Raw_Reads_Before", "Clean_Reads_After", "Total_Alignments", "Mapped_Reads"] if c in df_display.columns]
    for col in read_cols:
        df_display[col] = df_display[col].apply(format_reads)

    st.dataframe(
        df_display.style.apply(highlight_status, axis=1),
        use_container_width=True
    )

    if {"Clean_Reads_After", "Mapped_Reads", "Sample"}.issubset(df.columns):
        st.subheader("Reads Distribution")

        fig, ax = plt.subplots(figsize=(12, 5))
        plot_df = df.copy()
        x = range(len(plot_df))

        ax.bar(
            [i - 0.2 for i in x],
            plot_df["Clean_Reads_After"].apply(to_numeric_safe),
            width=0.4,
            label="Clean Reads"
        )
        ax.bar(
            [i + 0.2 for i in x],
            plot_df["Mapped_Reads"].apply(to_numeric_safe),
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

    low_samples = df[df["QC_Status"] == "Maybe Low Coverage"].copy()
    if len(low_samples) > 0:
        st.subheader("Flagged Samples")
        low_display = low_samples[["Sample", "Clean_Reads_After", "Mapped_Reads", "QC_Status"]].copy()
        if "Clean_Reads_After" in low_display.columns:
            low_display["Clean_Reads_After"] = low_display["Clean_Reads_After"].apply(format_reads)
        if "Mapped_Reads" in low_display.columns:
            low_display["Mapped_Reads"] = low_display["Mapped_Reads"].apply(format_reads)
        st.dataframe(low_display, use_container_width=True)
    else:
        st.success("All samples passed the minimum mapped reads threshold.")


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
# Tabs
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

    coverage_df = None

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

        df["QC_Status"] = df["Mapped_Reads"].apply(qc_status)
        df["Project"] = project_name_input

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

        run_info = {
            "Project": project_name_input,
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

        render_run_info(run_info)
        render_sample_qc(df)

        if coverage_file is not None:
            coverage_df = pd.read_csv(coverage_file, sep="\t")
            render_coverage_module(coverage_df)

# =========================
# SAVE CURRENT RUN
# =========================
st.subheader("Save to PostgreSQL")

if st.button("Save Current QC Run to PostgreSQL", type="primary"):
    try:
        run_record = normalize_record({
            "project_name": project_name_input,
            "flow_cell": flow_cell_input or "NA",
            "run_date_be": run_date_input or "NA",
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
        })

        sample_records = build_sample_records(df)
        coverage_records = build_coverage_records(coverage_df) if coverage_df is not None else []

        run_id = insert_qc_run(run_record, sample_records, coverage_records)
        st.success(f"QC run saved successfully. Run ID = {run_id}")

    except Exception as e:
        st.error(f"Failed to save QC run: {e}")

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
            pdf_df[col] = pdf_df[col].apply(format_reads)

        for col in ["% Duplication", "% > Q30", "GC content", "% Adapter", "Mapped_Reads(%)"]:
            pdf_df[col] = pdf_df[col].map(lambda x: f"{to_numeric_safe(x):.2f}")

        pdf_df["Mb Q30 bases"] = pdf_df["Mb Q30 bases"].map(lambda x: f"{to_numeric_safe(x):.3f}")
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


# =========================================================
# TAB 2: Browse Database
# =========================================================
with tab2:
    st.subheader("Project Selection")

    try:
        conn = get_connection()
        project_df = load_projects(conn)
    except Exception as e:
        st.error(f"Cannot connect to PostgreSQL database: {e}")
        st.stop()

    s1, s2 = st.columns([3, 1])
    with s1:
        project_search = st.text_input("Search Project", value="")
    with s2:
        search_clicked = st.button("Search")

    if search_clicked:
        st.session_state["project_search_value"] = project_search

    search_value = st.session_state.get("project_search_value", project_search)

    if search_value.strip():
        filtered_projects = project_df[
            project_df["project_name"].str.contains(search_value, case=False, na=False)
        ].copy()
    else:
        filtered_projects = project_df.copy()

    project_options = filtered_projects["project_name"].dropna().tolist()

    if len(project_options) == 0:
        st.warning("No matching project found.")
        conn.close()
        st.stop()

    selected_project = st.selectbox(
        "Select Project",
        project_options,
        key="db_selected_project"
    )

    run_df = load_runs_by_project(conn, selected_project)

    if run_df.empty:
        st.warning("No runs found for this project.")
        conn.close()
        st.stop()

    run_df["run_label"] = (
        run_df["flow_cell"].fillna("NA").astype(str) + " | " +
        run_df["run_date_be"].fillna("NA").astype(str) + " | Run ID: " +
        run_df["id"].astype(str)
    )

    selected_run_label = st.selectbox(
        "Select Run",
        run_df["run_label"].tolist(),
        key="db_selected_run"
    )

    selected_run_id = int(
        run_df.loc[run_df["run_label"] == selected_run_label, "id"].iloc[0]
    )

    run_info = load_run_info(conn, selected_run_id)
    sample_qc_df = load_sample_qc(conn, selected_run_id)
    coverage_df_db = load_target_coverage(conn, selected_run_id)

    conn.close()

    render_run_info(run_info)
    render_sample_qc(sample_qc_df)

    if not coverage_df_db.empty:
        render_coverage_module(coverage_df_db)
    
