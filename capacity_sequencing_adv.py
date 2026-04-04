import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import psycopg2
from psycopg2.extras import RealDictCursor
from datetime import datetime

st.set_page_config(page_title="NGS Capacity Planner", layout="wide")

# =========================
# CONSTANTS
# =========================
READ_LENGTH_BP = 300  # PE150 = 2 x 150
PARENTAGE_SIZE = 30645
WSSV_SIZE = 10603

KIT_OPTIONS = {
    "25M": {
        "capacity_bp": 7.5e9,
        "label": "25M (7.5 Gb)"
    },
    "60M": {
        "capacity_bp": 18e9,
        "label": "60M (18 Gb)"
    }
}

# =========================
# DB CONNECTION
# =========================
@st.cache_resource
def get_connection():
    return psycopg2.connect(
        host=st.secrets["postgres"]["host"],
        port=st.secrets["postgres"]["port"],
        dbname=st.secrets["postgres"]["dbname"],
        user=st.secrets["postgres"]["user"],
        password=st.secrets["postgres"]["password"],
        sslmode=st.secrets["postgres"].get("sslmode", "require")
    )

def insert_run(record: dict):
    conn = get_connection()
    with conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO sequencing_runs (
                    run_date,
                    reagent_kit,
                    reagent_kit_label,
                    coverage_x,
                    parentage_samples,
                    wssv_samples,
                    total_samples,
                    parentage_panel_bp,
                    wssv_panel_bp,
                    total_bp_required,
                    total_gb_required,
                    capacity_bp,
                    capacity_gb,
                    usage_percent,
                    remaining_bp,
                    remaining_gb,
                    reads_required,
                    reads_required_m,
                    reads_capacity,
                    reads_capacity_m,
                    reads_per_sample,
                    max_samples,
                    notes
                )
                VALUES (
                    %(run_date)s,
                    %(reagent_kit)s,
                    %(reagent_kit_label)s,
                    %(coverage_x)s,
                    %(parentage_samples)s,
                    %(wssv_samples)s,
                    %(total_samples)s,
                    %(parentage_panel_bp)s,
                    %(wssv_panel_bp)s,
                    %(total_bp_required)s,
                    %(total_gb_required)s,
                    %(capacity_bp)s,
                    %(capacity_gb)s,
                    %(usage_percent)s,
                    %(remaining_bp)s,
                    %(remaining_gb)s,
                    %(reads_required)s,
                    %(reads_required_m)s,
                    %(reads_capacity)s,
                    %(reads_capacity_m)s,
                    %(reads_per_sample)s,
                    %(max_samples)s,
                    %(notes)s
                )
                """,
                record
            )

def load_history(limit=100):
    conn = get_connection()
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(
            """
            SELECT
                id,
                run_date,
                reagent_kit,
                reagent_kit_label,
                coverage_x,
                parentage_samples,
                wssv_samples,
                total_samples,
                total_gb_required,
                capacity_gb,
                usage_percent,
                remaining_gb,
                reads_required_m,
                reads_capacity_m,
                reads_per_sample,
                max_samples,
                notes
            FROM sequencing_runs
            ORDER BY run_date DESC
            LIMIT %s
            """,
            (limit,)
        )
        rows = cur.fetchall()
    return pd.DataFrame(rows)

# =========================
# HELPERS
# =========================
def calculate_metrics(kit_key, coverage, parentage_samples, wssv_samples):
    capacity_bp = KIT_OPTIONS[kit_key]["capacity_bp"]
    reagent_kit_label = KIT_OPTIONS[kit_key]["label"]

    parentage_bp = parentage_samples * PARENTAGE_SIZE * coverage
    wssv_bp = wssv_samples * WSSV_SIZE * coverage
    total_bp = parentage_bp + wssv_bp

    total_gb = total_bp / 1e9
    capacity_gb = capacity_bp / 1e9

    reads_required = total_bp / READ_LENGTH_BP
    reads_capacity = capacity_bp / READ_LENGTH_BP

    usage_percent = (total_bp / capacity_bp * 100) if capacity_bp > 0 else 0
    remaining_bp = capacity_bp - total_bp
    remaining_gb = remaining_bp / 1e9

    total_samples = parentage_samples + wssv_samples
    reads_per_sample = (reads_required / total_samples) if total_samples > 0 else 0

    avg_bp_per_sample = (total_bp / total_samples) if total_samples > 0 else 0
    max_samples = (capacity_bp / avg_bp_per_sample) if avg_bp_per_sample > 0 else 0

    return {
        "reagent_kit": kit_key,
        "reagent_kit_label": reagent_kit_label,
        "coverage_x": coverage,
        "parentage_samples": parentage_samples,
        "wssv_samples": wssv_samples,
        "total_samples": total_samples,
        "parentage_panel_bp": PARENTAGE_SIZE,
        "wssv_panel_bp": WSSV_SIZE,
        "parentage_bp_required": parentage_bp,
        "wssv_bp_required": wssv_bp,
        "total_bp_required": total_bp,
        "total_gb_required": total_gb,
        "capacity_bp": capacity_bp,
        "capacity_gb": capacity_gb,
        "usage_percent": usage_percent,
        "remaining_bp": remaining_bp,
        "remaining_gb": remaining_gb,
        "reads_required": reads_required,
        "reads_required_m": reads_required / 1e6,
        "reads_capacity": reads_capacity,
        "reads_capacity_m": reads_capacity / 1e6,
        "reads_per_sample": reads_per_sample,
        "max_samples": max_samples,
    }

# =========================
# TITLE
# =========================
st.title("🧬 NGS Capacity Planner")
st.caption("Salus Nimbo planning tool with PostgreSQL run history")

# =========================
# LOAD HISTORY FIRST
# =========================
st.subheader("📚 Historical Runs")
history_df = load_history(limit=200)

if not history_df.empty:
    display_df = history_df.copy()
    if "run_date" in display_df.columns:
        display_df["run_date"] = pd.to_datetime(display_df["run_date"]).dt.strftime("%m/%d/%Y %H:%M")

    st.dataframe(display_df, use_container_width=True)

    st.download_button(
        "Download History CSV",
        history_df.to_csv(index=False),
        file_name="sequencing_run_history.csv",
        mime="text/csv"
    )
else:
    st.info("No saved run history yet.")

st.divider()

# =========================
# INPUTS
# =========================
st.subheader("⚙️ Current Run Setup")

col1, col2, col3 = st.columns(3)

with col1:
    kit_key = st.radio(
        "Select Reagent Kit",
        options=list(KIT_OPTIONS.keys()),
        format_func=lambda x: KIT_OPTIONS[x]["label"],
        horizontal=True
    )

with col2:
    coverage = st.number_input("Coverage (X)", min_value=1, max_value=1000, value=20, step=1)

with col3:
    notes = st.text_input("Notes / Project name", value="")

col4, col5 = st.columns(2)

with col4:
    parentage_samples = st.number_input("Parentage samples", min_value=0, max_value=100000, value=300, step=1)

with col5:
    wssv_samples = st.number_input("WSSV samples", min_value=0, max_value=100000, value=200, step=1)

# =========================
# CALCULATE
# =========================
metrics = calculate_metrics(
    kit_key=kit_key,
    coverage=coverage,
    parentage_samples=parentage_samples,
    wssv_samples=wssv_samples
)

# =========================
# SUMMARY METRICS
# =========================
st.subheader("📊 Current Run Summary")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Total Required (Gb)", f"{metrics['total_gb_required']:.4f}")
m2.metric("Capacity (Gb)", f"{metrics['capacity_gb']:.2f}")
m3.metric("Usage (%)", f"{metrics['usage_percent']:.2f}")
m4.metric("Remaining (Gb)", f"{metrics['remaining_gb']:.4f}")

m5, m6, m7, m8 = st.columns(4)
m5.metric("Reads Required (M)", f"{metrics['reads_required_m']:.4f}")
m6.metric("Reads Capacity (M)", f"{metrics['reads_capacity_m']:.2f}")
m7.metric("Reads / Sample", f"{metrics['reads_per_sample']:.2f}")
m8.metric("Max Samples / Run", f"{metrics['max_samples']:.0f}")

if metrics["total_bp_required"] > metrics["capacity_bp"]:
    st.error("⚠️ This setup exceeds sequencing capacity.")
else:
    st.success("✅ This setup is within sequencing capacity.")

# =========================
# PLOTS
# =========================
plot_col1, plot_col2 = st.columns(2)

with plot_col1:
    fig_gb = go.Figure()
    fig_gb.add_trace(go.Bar(
        x=["Parentage", "WSSV"],
        y=[
            metrics["parentage_bp_required"] / 1e9,
            metrics["wssv_bp_required"] / 1e9
        ],
        name="Required Gb"
    ))
    fig_gb.add_hline(y=metrics["capacity_gb"])
    fig_gb.update_layout(
        title="Sequencing Usage (Gb)",
        yaxis_title="Gb"
    )
    st.plotly_chart(fig_gb, use_container_width=True)

with plot_col2:
    fig_reads = go.Figure()
    fig_reads.add_trace(go.Bar(
        x=["Required Reads"],
        y=[metrics["reads_required_m"]],
        name="Required Reads (M)"
    ))
    fig_reads.add_trace(go.Bar(
        x=["Capacity Reads"],
        y=[metrics["reads_capacity_m"]],
        name="Capacity Reads (M)"
    ))
    fig_reads.update_layout(
        title="Reads Usage (Million Reads)",
        yaxis_title="Reads (M)"
    )
    st.plotly_chart(fig_reads, use_container_width=True)

# =========================
# CURRENT RUN EXPORT
# =========================
st.subheader("📁 Export Current Run")

current_run_df = pd.DataFrame([{
    "Run_Date": datetime.now().strftime("%m/%d/%Y %H:%M"),
    "Selected_Reagent_Kit": metrics["reagent_kit"],
    "Selected_Reagent_Kit_Label": metrics["reagent_kit_label"],
    "Coverage_X": metrics["coverage_x"],
    "Parentage_Samples": metrics["parentage_samples"],
    "WSSV_Samples": metrics["wssv_samples"],
    "Total_Samples": metrics["total_samples"],
    "Parentage_Panel_bp": metrics["parentage_panel_bp"],
    "WSSV_Panel_bp": metrics["wssv_panel_bp"],
    "Total_Gb_Required": metrics["total_gb_required"],
    "Capacity_Gb": metrics["capacity_gb"],
    "Usage_Percent": metrics["usage_percent"],
    "Remaining_Gb": metrics["remaining_gb"],
    "Reads_Required_M": metrics["reads_required_m"],
    "Reads_Capacity_M": metrics["reads_capacity_m"],
    "Reads_per_Sample": metrics["reads_per_sample"],
    "Max_Samples": metrics["max_samples"],
    "Notes": notes
}])

st.dataframe(current_run_df, use_container_width=True)

st.download_button(
    "Download Current Run CSV",
    current_run_df.to_csv(index=False),
    file_name="NGS_planning_report.csv",
    mime="text/csv"
)

# =========================
# SAVE TO DATABASE
# =========================
st.subheader("💾 Save Current Run")

if st.button("Save Current Run to PostgreSQL", type="primary"):
    record = {
        "run_date": datetime.now(),
        "reagent_kit": metrics["reagent_kit"],
        "reagent_kit_label": metrics["reagent_kit_label"],
        "coverage_x": metrics["coverage_x"],
        "parentage_samples": metrics["parentage_samples"],
        "wssv_samples": metrics["wssv_samples"],
        "total_samples": metrics["total_samples"],
        "parentage_panel_bp": metrics["parentage_panel_bp"],
        "wssv_panel_bp": metrics["wssv_panel_bp"],
        "total_bp_required": metrics["total_bp_required"],
        "total_gb_required": metrics["total_gb_required"],
        "capacity_bp": metrics["capacity_bp"],
        "capacity_gb": metrics["capacity_gb"],
        "usage_percent": metrics["usage_percent"],
        "remaining_bp": metrics["remaining_bp"],
        "remaining_gb": metrics["remaining_gb"],
        "reads_required": metrics["reads_required"],
        "reads_required_m": metrics["reads_required_m"],
        "reads_capacity": metrics["reads_capacity"],
        "reads_capacity_m": metrics["reads_capacity_m"],
        "reads_per_sample": metrics["reads_per_sample"],
        "max_samples": metrics["max_samples"],
        "notes": notes
    }

    try:
        insert_run(record)
        st.success("Run saved successfully.")
        st.rerun()
    except Exception as e:
        st.error(f"Failed to save run: {e}")
