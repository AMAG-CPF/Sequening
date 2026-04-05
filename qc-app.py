import io
import re
import hmac
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2
import streamlit as st
from matplotlib.ticker import FuncFormatter
from psycopg2.extras import RealDictCursor


# =========================================================
# Streamlit config
# =========================================================
st.set_page_config(page_title="NGS QC Dashboard", layout="wide")
st.title("🧬 NGS QC Dashboard")


# =========================================================
# SECURITY
# =========================================================
def check_password():
    def password_entered():
        expected = st.secrets["app"]["password"]
        entered = st.session_state.get("password", "")
        st.session_state["password_correct"] = hmac.compare_digest(entered, expected)

    if st.session_state.get("password_correct", False):
        return

    st.text_input("App Password", type="password", key="password", on_change=password_entered)

    if st.session_state.get("password"):
        st.error("Incorrect password")

    st.stop()


check_password()


# =========================================================
# DATABASE
# =========================================================
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
    elif report_df.s
