import streamlit as st
import plotly.graph_objects as go
from protein_fetcher import fetch_pdb, fetch_protein_info
from visualizer import render_protein

# ── Page config ─────────────────────────────────────
st.set_page_config(
    page_title="Protein Explorer",
    page_icon="🔬",
    layout="wide"
)

# ── Header ───────────────────────────────────────────
st.title("🔬 Protein Structure Explorer")
st.markdown("Visualize any protein in 3D using AlphaFold's database")

# ── Sidebar — controls ───────────────────────────────
with st.sidebar:
    st.header("🧬 Load a Protein")

    # Quick presets
    st.markdown("**Try these examples:**")
    examples = {
        "Hemoglobin (Sickle Cell)": "P69905",
        "Insulin":                  "P01308",
        "BRCA1 (Breast Cancer)":    "P38398",
        "Spike Protein (COVID-19)": "P0DTC2",
        "Collagen Type I":          "P02452",
    }

    for label, uid in examples.items():
        if st.button(label, use_container_width=True):
            st.session_state["uniprot_id"] = uid

    st.divider()

    # Manual input
    uniprot_id = st.text_input(
        "Or enter UniProt ID manually:",
        value=st.session_state.get("uniprot_id", "P69905")
    ).strip().upper()

    load_btn = st.button("🔍 Load Protein", type="primary", use_container_width=True)

    st.divider()

    # Display options
    st.header("🎨 Display Options")
    style = st.selectbox("Representation", ["cartoon", "stick", "sphere", "surface"])
    color = st.selectbox("Color Scheme",   ["confidence", "rainbow", "secondary"])

# ── Main content ─────────────────────────────────────
if load_btn or "uniprot_id" in st.session_state:

    uid = st.session_state.get("uniprot_id", uniprot_id)
    if load_btn:
        uid = uniprot_id

    with st.spinner(f"Fetching protein {uid} from AlphaFold..."):
        pdb_data = fetch_pdb(uid)
        info     = fetch_protein_info(uid)

    if pdb_data is None:
        st.error(f"❌ Could not find protein '{uid}'. Check the UniProt ID.")
    else:
        # ── Two-column layout ─────────────────────
        col1, col2 = st.columns([2, 1])

        with col1:
            st.subheader("🧬 3D Structure")
            render_protein(pdb_data, style=style, color=color)

            # Color legend for confidence view
            if color == "confidence":
                st.markdown("""
                **pLDDT Confidence Score:**
                🔵 **Blue (90–100)** = Very high confidence &nbsp;|&nbsp;
                🟢 **Green (70–90)** = Confident &nbsp;|&nbsp;
                🟡 **Yellow (50–70)** = Low confidence &nbsp;|&nbsp;
                🔴 **Red (<50)** = Very low confidence
                """)

        with col2:
            if info:
                st.subheader("📋 Protein Info")
                st.metric("Protein Name", info["name"])
                st.metric("Organism",     info["organism"])
                st.metric("Length",       f"{info['length']} amino acids")

                st.markdown("**Function:**")
                st.info(info["function"][:400] + "..." if len(info["function"]) > 400 else info["function"])

                # External links
                st.markdown("**🔗 Learn More:**")
                st.markdown(f"[UniProt Page](https://www.uniprot.org/uniprot/{uid})")
                st.markdown(f"[AlphaFold Page](https://alphafold.ebi.ac.uk/entry/{uid})")

        # ── Confidence score chart ────────────────
        st.subheader("📊 Confidence Score Per Residue (pLDDT)")
        st.markdown("Shows how confident AlphaFold is about each part of the structure")

        # Parse B-factor (pLDDT) from PDB
        residues, scores = [], []
        seen = set()
        for line in pdb_data.split("\n"):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_num = int(line[22:26].strip())
                if res_num not in seen:
                    try:
                        b_factor = float(line[60:66].strip())
                        residues.append(res_num)
                        scores.append(b_factor)
                        seen.add(res_num)
                    except:
                        pass

        if residues:
            colors = ["#003f7f" if s >= 90 else "#1f9e89" if s >= 70 else "#f5a623" if s >= 50 else "#d73027"
                      for s in scores]

            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=residues, y=scores,
                marker_color=colors,
                hovertemplate="Residue %{x}<br>pLDDT: %{y:.1f}<extra></extra>"
            ))
            fig.add_hline(y=90, line_dash="dot", line_color="#003f7f", annotation_text="Very High (90)")
            fig.add_hline(y=70, line_dash="dot", line_color="#1f9e89", annotation_text="Confident (70)")
            fig.add_hline(y=50, line_dash="dot", line_color="#f5a623", annotation_text="Low (50)")

            fig.update_layout(
                height=300,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis_title="Residue Position",
                yaxis_title="pLDDT Score",
                yaxis_range=[0, 100],
                showlegend=False
            )
            st.plotly_chart(fig, use_container_width=True)

else:
    # Landing state — nothing loaded yet
    st.info("👈 Select a protein from the sidebar or enter a UniProt ID to get started!")

    st.markdown("### 🧪 What is a UniProt ID?")
    st.markdown("""
    Every protein has a unique ID in the UniProt database — like a passport number.
    For example:
    - `P69905` → Human Hemoglobin (carries oxygen in blood)
    - `P01308` → Human Insulin (regulates blood sugar)
    - `P38398` → BRCA1 (breast cancer related gene)

    Just enter any ID and the app fetches the 3D structure automatically!
    """)