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
READ_LENGTH_BP = 300  # PE150 = 2x150 bp

KIT_OPTIONS = {
    "25M": {"capacity_bp": 7.5e9, "label": "25M (7.5 Gb)"},
    "60M": {"capacity_bp": 18e9, "label": "60M (18 Gb)"}
}

PRESET_PANELS = {
    "Parentage": 30645,
    "WSSV": 10603
}

# =========================
# DB
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

def insert_run(summary_record: dict, panel_records: list[dict]) -> int:
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
                    total_panels,
                    total_samples,
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
                    %(total_panels)s,
                    %(total_samples)s,
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
                RETURNING id
                """,
                summary_record
            )
            run_id = cur.fetchone()[0]

            for panel in panel_records:
                panel["run_id"] = run_id
                cur.execute(
                    """
                    INSERT INTO sequencing_run_panels (
                        run_id,
                        panel_name,
                        panel_size_bp,
                        samples,
                        coverage_x,
                        required_bp,
                        required_gb
                    )
                    VALUES (
                        %(run_id)s,
                        %(panel_name)s,
                        %(panel_size_bp)s,
                        %(samples)s,
                        %(coverage_x)s,
                        %(required_bp)s,
                        %(required_gb)s
                    )
                    """,
                    panel
                )
    return run_id

def load_run_history(limit=100):
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
                total_panels,
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

def load_panel_history(limit_runs=20):
    conn = get_connection()
    with conn.cursor(cursor_factory=RealDictCursor) as cur:
        cur.execute(
            """
            SELECT
                r.id AS run_id,
                r.run_date,
                r.reagent_kit,
                r.coverage_x,
                p.panel_name,
                p.panel_size_bp,
                p.samples,
                p.required_gb
            FROM sequencing_runs r
            JOIN sequencing_run_panels p
              ON r.id = p.run_id
            ORDER BY r.run_date DESC, p.id ASC
            LIMIT %s
            """,
            (limit_runs * 20,)
        )
        rows = cur.fetchall()
    return pd.DataFrame(rows)

# =========================
# SESSION STATE
# =========================
if "custom_panels" not in st.session_state:
    st.session_state.custom_panels = []

def add_custom_panel():
    st.session_state.custom_panels.append({
        "name": "",
        "size_bp": 1,
        "samples": 0
    })

def remove_custom_panel(index: int):
    st.session_state.custom_panels.pop(index)

# =========================
# HELPERS
# =========================
def build_selected_panels(coverage: int):
    selected_panels = []

    for panel_name, default_size in PRESET_PANELS.items():
        use_panel = st.session_state.get(f"use_{panel_name}", False)
        if use_panel:
            samples = st.session_state.get(f"{panel_name}_samples", 0)
            if samples > 0:
                selected_panels.append({
                    "panel_name": panel_name,
                    "panel_size_bp": default_size,
                    "samples": samples,
                    "coverage_x": coverage
                })

    for panel in st.session_state.custom_panels:
        name = str(panel["name"]).strip()
        size_bp = int(panel["size_bp"])
        samples = int(panel["samples"])
        if name and size_bp > 0 and samples > 0:
            selected_panels.append({
                "panel_name": name,
                "panel_size_bp": size_bp,
                "samples": samples,
                "coverage_x": coverage
            })

    return selected_panels

def calculate_metrics(kit_key: str, coverage: int, panels: list[dict]):
    capacity_bp = KIT_OPTIONS[kit_key]["capacity_bp"]
    capacity_gb = capacity_bp / 1e9

    panel_rows = []
    total_bp_required = 0
    total_samples = 0

    for panel in panels:
        required_bp = panel["panel_size_bp"] * panel["samples"] * coverage
        required_gb = required_bp / 1e9

        total_bp_required += required_bp
        total_samples += panel["samples"]

        panel_rows.append({
            "Panel": panel["panel_name"],
            "Panel_Size_bp": panel["panel_size_bp"],
            "Samples": panel["samples"],
            "Coverage_X": coverage,
            "Required_bp": required_bp,
            "Required_Gb": required_gb
        })

    total_gb_required = total_bp_required / 1e9
    reads_required = total_bp_required / READ_LENGTH_BP
    reads_capacity = capacity_bp / READ_LENGTH_BP

    usage_percent = (total_bp_required / capacity_bp * 100) if capacity_bp > 0 else 0
    remaining_bp = capacity_bp - total_bp_required
    remaining_gb = remaining_bp / 1e9
    reads_per_sample = (reads_required / total_samples) if total_samples > 0 else 0

    avg_bp_per_sample = (total_bp_required / total_samples) if total_samples > 0 else 0
    max_samples = (capacity_bp / avg_bp_per_sample) if avg_bp_per_sample > 0 else 0

    panels_df = pd.DataFrame(panel_rows)

    return {
        "reagent_kit": kit_key,
        "reagent_kit_label": KIT_OPTIONS[kit_key]["label"],
        "coverage_x": coverage,
        "total_panels": len(panels),
        "total_samples": total_samples,
        "total_bp_required": total_bp_required,
        "total_gb_required": total_gb_required,
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
        "panels_df": panels_df
    }

# =========================
# TITLE
# =========================
st.title("🧬 NGS Capacity Planner")
st.caption("Dynamic panel planning with PostgreSQL run history")

# =========================
# HISTORY
# =========================
st.subheader("📚 Historical Runs")
history_df = load_run_history(limit=200)

if not history_df.empty:
    hist_display = history_df.copy()
    hist_display["run_date"] = pd.to_datetime(hist_display["run_date"]).dt.strftime("%m/%d/%Y %H:%M")
    st.dataframe(hist_display, use_container_width=True)

    st.download_button(
        "Download Run History CSV",
        history_df.to_csv(index=False),
        file_name="sequencing_run_history.csv",
        mime="text/csv"
    )
else:
    st.info("No saved run history yet.")

with st.expander("Show historical panel breakdown"):
    panel_history_df = load_panel_history(limit_runs=50)
    if not panel_history_df.empty:
        panel_history_display = panel_history_df.copy()
        panel_history_display["run_date"] = pd.to_datetime(panel_history_display["run_date"]).dt.strftime("%m/%d/%Y %H:%M")
        st.dataframe(panel_history_display, use_container_width=True)

        st.download_button(
            "Download Panel History CSV",
            panel_history_df.to_csv(index=False),
            file_name="sequencing_panel_history.csv",
            mime="text/csv"
        )
    else:
        st.info("No saved panel history yet.")

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

st.markdown("### Preset Panels")

for panel_name, panel_size in PRESET_PANELS.items():
    c1, c2, c3 = st.columns([1.4, 1.2, 1.2])
    with c1:
        st.checkbox(
            f"Use {panel_name}",
            value=(panel_name == "Parentage"),
            key=f"use_{panel_name}"
        )
    with c2:
        st.number_input(
            f"{panel_name} size (bp)",
            min_value=1,
            value=panel_size,
            step=1,
            disabled=True,
            key=f"{panel_name}_size_display"
        )
    with c3:
        st.number_input(
            f"{panel_name} samples",
            min_value=0,
            max_value=100000,
            value=300 if panel_name == "Parentage" else 0,
            step=1,
            key=f"{panel_name}_samples"
        )

st.markdown("### Custom Panels")

if st.button("Add Custom Panel"):
    add_custom_panel()

for i, panel in enumerate(st.session_state.custom_panels):
    st.markdown(f"**Custom Panel {i+1}**")
    c1, c2, c3, c4 = st.columns([1.5, 1.2, 1.1, 0.8])

    with c1:
        panel["name"] = st.text_input(
            f"Panel name #{i+1}",
            value=panel["name"],
            key=f"custom_name_{i}"
        )

    with c2:
        panel["size_bp"] = st.number_input(
            f"Panel size (bp) #{i+1}",
            min_value=1,
            value=int(panel["size_bp"]),
            step=1,
            key=f"custom_size_{i}"
        )

    with c3:
        panel["samples"] = st.number_input(
            f"Samples #{i+1}",
            min_value=0,
            value=int(panel["samples"]),
            step=1,
            key=f"custom_samples_{i}"
        )

    with c4:
        st.write("")
        st.write("")
        if st.button("Remove", key=f"remove_panel_{i}"):
            remove_custom_panel(i)
            st.rerun()

selected_panels = build_selected_panels(coverage=coverage)
metrics = calculate_metrics(kit_key=kit_key, coverage=coverage, panels=selected_panels)

# =========================
# SUMMARY
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
# PANEL BREAKDOWN
# =========================
st.subheader("🧾 Panel Breakdown")

if not metrics["panels_df"].empty:
    st.dataframe(metrics["panels_df"], use_container_width=True)
else:
    st.info("No valid panel selected yet.")

# =========================
# PLOTS
# =========================
if not metrics["panels_df"].empty:
    colp1, colp2 = st.columns(2)

    with colp1:
        fig_gb = go.Figure()
        fig_gb.add_trace(go.Bar(
            x=metrics["panels_df"]["Panel"],
            y=metrics["panels_df"]["Required_Gb"],
            name="Required Gb"
        ))
        fig_gb.add_hline(y=metrics["capacity_gb"])
        fig_gb.update_layout(
            title="Sequencing Usage by Panel (Gb)",
            yaxis_title="Gb"
        )
        st.plotly_chart(fig_gb, use_container_width=True)

    with colp2:
        fig_reads = go.Figure()
        fig_reads.add_trace(go.Bar(
            x=["Required Reads", "Capacity Reads"],
            y=[metrics["reads_required_m"], metrics["reads_capacity_m"]],
            name="Reads (M)"
        ))
        fig_reads.update_layout(
            title="Reads Usage (Million Reads)",
            yaxis_title="Reads (M)"
        )
        st.plotly_chart(fig_reads, use_container_width=True)

# =========================
# EXPORT
# =========================
st.subheader("📁 Export Current Run")

summary_df = pd.DataFrame([{
    "Run_Date": datetime.now().strftime("%m/%d/%Y %H:%M"),
    "Selected_Reagent_Kit": metrics["reagent_kit"],
    "Selected_Reagent_Kit_Label": metrics["reagent_kit_label"],
    "Coverage_X": metrics["coverage_x"],
    "Total_Panels": metrics["total_panels"],
    "Total_Samples": metrics["total_samples"],
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

st.dataframe(summary_df, use_container_width=True)

st.download_button(
    "Download Summary CSV",
    summary_df.to_csv(index=False),
    file_name="NGS_planning_summary.csv",
    mime="text/csv"
)

if not metrics["panels_df"].empty:
    st.download_button(
        "Download Panel Breakdown CSV",
        metrics["panels_df"].to_csv(index=False),
        file_name="NGS_panel_breakdown.csv",
        mime="text/csv"
    )

# =========================
# SAVE
# =========================
st.subheader("💾 Save Current Run")

if st.button("Save Current Run to PostgreSQL", type="primary"):
    if metrics["panels_df"].empty:
        st.warning("Please select at least one valid panel with samples > 0 before saving.")
    else:
        summary_record = {
            "run_date": datetime.now(),
            "reagent_kit": metrics["reagent_kit"],
            "reagent_kit_label": metrics["reagent_kit_label"],
            "coverage_x": metrics["coverage_x"],
            "total_panels": metrics["total_panels"],
            "total_samples": metrics["total_samples"],
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

        panel_records = []
        for _, row in metrics["panels_df"].iterrows():
            panel_records.append({
                "panel_name": row["Panel"],
                "panel_size_bp": int(row["Panel_Size_bp"]),
                "samples": int(row["Samples"]),
                "coverage_x": int(row["Coverage_X"]),
                "required_bp": float(row["Required_bp"]),
                "required_gb": float(row["Required_Gb"])
            })

        try:
            run_id = insert_run(summary_record, panel_records)
            st.success(f"Run saved successfully. Run ID = {run_id}")
            st.rerun()
        except Exception as e:
            st.error(f"Failed to save run: {e}")
