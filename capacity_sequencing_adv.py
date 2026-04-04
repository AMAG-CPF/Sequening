import io
import re
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from matplotlib.ticker import FuncFormatter
import psycopg2
import pandas as pd
from psycopg2.extras import RealDictCursor

@st.cache_resource
def get_connection():
    return psycopg2.connect(
        host="localhost",
        port=5432,
        dbname="sequencing_app",
        user="postgres",
        password="AMAG"
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

def load_panel_history(limit=500):
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
            (limit,)
        )
        rows = cur.fetchall()
    return pd.DataFrame(rows)


st.set_page_config(page_title="NGS Capacity Planner", layout="wide")
st.subheader("📚 Historical Runs")
READ_LENGTH_BP = 300  # PE150

KIT_OPTIONS = {
    "25M": {"capacity_bp": 7.5e9, "label": "25M (7.5 Gb)"},
    "60M": {"capacity_bp": 18e9, "label": "60M (18 Gb)"}
}

PRESET_PANELS = {
    "Parentage": 30645,
    "WSSV": 10603
}

# -------------------------
# session state for custom panels
# -------------------------
if "custom_panels" not in st.session_state:
    st.session_state.custom_panels = []

def add_custom_panel():
    st.session_state.custom_panels.append({
        "name": "",
        "size_bp": 0,
        "samples": 0
    })
def remove_custom_panel(index: int):
    st.session_state.custom_panels.pop(index)
st.title(":material/home: NGS Capacity Planner")

# -------------------------
# HISTORY placeholder
# -------------------------
st.subheader("⚙️ Current Run Setup")

col1, col2 = st.columns(2)
with col1:
    kit_key = st.radio(
        "Select Reagent Kit",
        options=list(KIT_OPTIONS.keys()),
        format_func=lambda x: KIT_OPTIONS[x]["label"],
        horizontal=True
    )
with col2:
    coverage = st.number_input("Coverage (X)", min_value=1, max_value=1000, value=20, step=1)

notes = st.text_input("Notes / Project name", value="")

st.markdown("### Preset Panels")

selected_panels = []

for panel_name, panel_size in PRESET_PANELS.items():
    use_panel = st.checkbox(f"Use {panel_name} ({panel_size:,} bp)", value=(panel_name in ["Parentage", "WSSV"]))
    if use_panel:
        c1, c2 = st.columns([2, 1])
        with c1:
            samples = st.number_input(
                f"{panel_name} samples",
                min_value=0,
                max_value=100000,
                value=300 if panel_name == "Parentage" else 200,
                step=1,
                key=f"{panel_name}_samples"
            )
        with c2:
            st.number_input(
                f"{panel_name} panel size (bp)",
                min_value=1,
                value=panel_size,
                step=1,
                key=f"{panel_name}_size",
                disabled=True
            )

        selected_panels.append({
            "panel_name": panel_name,
            "panel_size_bp": panel_size,
            "samples": samples
        })

st.markdown("### Custom Panels")

c_add, _ = st.columns([1, 4])
with c_add:
    if st.button("Add Custom Panel"):
        add_custom_panel()

for i, panel in enumerate(st.session_state.custom_panels):
    st.markdown(f"**Custom Panel {i+1}**")
    c1, c2, c3, c4 = st.columns(4)

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
            value=max(1, int(panel["size_bp"])) if panel["size_bp"] else 1,
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

    if panel["name"].strip():
        selected_panels.append({
            "panel_name": panel["name"].strip(),
            "panel_size_bp": panel["size_bp"],
            "samples": panel["samples"]
        })
    with c4:
        st.write("")  # spacing
        st.write("")
        if st.button(f"Remove Panel #{i+1}", key=f"remove_{i}"):
            st.session_state.custom_panels.pop(i)
            st.rerun()
# -------------------------
# calculations
# -------------------------
capacity_bp = KIT_OPTIONS[kit_key]["capacity_bp"]
capacity_gb = capacity_bp / 1e9

panel_rows = []
total_bp_required = 0
total_samples = 0

for panel in selected_panels:
    panel_bp_required = panel["panel_size_bp"] * panel["samples"] * coverage
    total_bp_required += panel_bp_required
    total_samples += panel["samples"]

    panel_rows.append({
        "Panel": panel["panel_name"],
        "Panel_Size_bp": panel["panel_size_bp"],
        "Samples": panel["samples"],
        "Coverage_X": coverage,
        "Required_bp": panel_bp_required,
        "Required_Gb": panel_bp_required / 1e9
    })

panels_df = pd.DataFrame(panel_rows)

reads_required = total_bp_required / READ_LENGTH_BP
reads_capacity = capacity_bp / READ_LENGTH_BP
usage_percent = (total_bp_required / capacity_bp * 100) if capacity_bp > 0 else 0
remaining_bp = capacity_bp - total_bp_required
remaining_gb = remaining_bp / 1e9
reads_per_sample = (reads_required / total_samples) if total_samples > 0 else 0
avg_bp_per_sample = (total_bp_required / total_samples) if total_samples > 0 else 0
max_samples = (capacity_bp / avg_bp_per_sample) if avg_bp_per_sample > 0 else 0

# -------------------------
# summary
# -------------------------
st.subheader("📊 Current Run Summary")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Total Required (Gb)", f"{total_bp_required/1e9:.4f}")
m2.metric("Capacity (Gb)", f"{capacity_gb:.2f}")
m3.metric("Usage (%)", f"{usage_percent:.2f}")
m4.metric("Remaining (Gb)", f"{remaining_gb:.4f}")

m5, m6, m7, m8 = st.columns(4)
m5.metric("Reads Required (M)", f"{reads_required/1e6:.4f}")
m6.metric("Reads Capacity (M)", f"{reads_capacity/1e6:.2f}")
m7.metric("Reads / Sample", f"{reads_per_sample:.2f}")
m8.metric("Max Samples / Run", f"{max_samples:.0f}")

if total_bp_required > capacity_bp:
    st.error("⚠️ This setup exceeds sequencing capacity.")
else:
    st.success("✅ This setup is within sequencing capacity.")

# -------------------------
# panel table
# -------------------------
st.subheader("🧾 Panel Breakdown")
if not panels_df.empty:
    st.dataframe(panels_df, use_container_width=True)
else:
    st.info("No panels selected yet.")

# -------------------------
# plots
# -------------------------
if not panels_df.empty:
    fig_gb = go.Figure()
    fig_gb.add_trace(go.Bar(
        x=panels_df["Panel"],
        y=panels_df["Required_Gb"],
        name="Required Gb"
    ))
    fig_gb.add_hline(y=capacity_gb)
    fig_gb.update_layout(
        title="Sequencing Usage by Panel (Gb)",
        yaxis_title="Gb"
    )
    st.plotly_chart(fig_gb, use_container_width=True)

# -------------------------
# current run export
# -------------------------
st.subheader("📁 Export Current Run")

summary_df = pd.DataFrame([{
    "Run_Date": datetime.now().strftime("%m/%d/%Y %H:%M"),
    "Selected_Reagent_Kit": kit_key,
    "Selected_Reagent_Kit_Label": KIT_OPTIONS[kit_key]["label"],
    "Coverage_X": coverage,
    "Total_Panels": len(selected_panels),
    "Total_Samples": total_samples,
    "Total_Gb_Required": total_bp_required / 1e9,
    "Capacity_Gb": capacity_gb,
    "Usage_Percent": usage_percent,
    "Remaining_Gb": remaining_gb,
    "Reads_Required_M": reads_required / 1e6,
    "Reads_Capacity_M": reads_capacity / 1e6,
    "Reads_per_Sample": reads_per_sample,
    "Max_Samples": max_samples,
    "Notes": notes
}])

st.dataframe(summary_df, use_container_width=True)

csv_detail = panels_df.to_csv(index=False)
csv_summary = summary_df.to_csv(index=False)

st.download_button(
    "Download Summary CSV",
    csv_summary,
    file_name="NGS_planning_summary.csv",
    mime="text/csv"
)

st.download_button(
    "Download Panel Breakdown CSV",
    csv_detail,
    file_name="NGS_panel_breakdown.csv",
    mime="text/csv"
)

st.subheader("💾 Save Current Run")

if st.button("Save Current Run to PostgreSQL", type="primary"):
    if panels_df.empty:
        st.warning("Please add at least one valid panel before saving.")
    else:
        summary_record = {
            "run_date": datetime.now(),
            "reagent_kit": kit_key,
            "reagent_kit_label": KIT_OPTIONS[kit_key]["label"],
            "coverage_x": coverage,
            "total_panels": len(selected_panels),
            "total_samples": total_samples,
            "total_bp_required": total_bp_required,
            "total_gb_required": total_bp_required / 1e9,
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
            "notes": notes
        }

        panel_records = []
        for _, row in panels_df.iterrows():
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
