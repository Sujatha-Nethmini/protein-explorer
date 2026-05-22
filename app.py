import streamlit as st
import plotly.graph_objects as go
from protein_fetcher import fetch_pdb, fetch_protein_info
from visualizer import render_protein, render_protein_with_mutation
from mutation_analyzer import (
    get_residue_at_position,
    get_plddt_at_position,
    calculate_ddg,
    interpret_ddg
)

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(page_title="Protein Explorer", page_icon="🔬", layout="wide")

# ── Initialize session state ─────────────────────────────────────────────────
# This is the fix — store everything here so it survives widget interactions
if "pdb_data" not in st.session_state:
    st.session_state.pdb_data = None
if "protein_info" not in st.session_state:
    st.session_state.protein_info = None
if "loaded_id" not in st.session_state:
    st.session_state.loaded_id = ""
if "mutation_result" not in st.session_state:
    st.session_state.mutation_result = None

# ── Header ───────────────────────────────────────────────────────────────────
st.title("🔬 Protein Structure Explorer")
st.markdown("Visualize any protein in 3D using AlphaFold's database")

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🧬 Load a Protein")

    examples = {
        "Hemoglobin (Sickle Cell)": "P69905",
        "Insulin":                  "P01308",
        "BRCA1 (Breast Cancer)":    "P38398",
        "Spike Protein (COVID-19)": "P0DTC2",
        "Collagen Type I":          "P02452",
    }

    st.markdown("**Try these examples:**")
    for label, uid in examples.items():
        if st.button(label, use_container_width=True):
            st.session_state.loaded_id = uid
            st.session_state.pdb_data = None       # force reload
            st.session_state.protein_info = None
            st.session_state.mutation_result = None

    st.divider()

    uniprot_input = st.text_input(
        "Or enter UniProt ID manually:",
        value=st.session_state.loaded_id
    ).strip().upper()

    if st.button("🔍 Load Protein", type="primary", use_container_width=True):
        st.session_state.loaded_id = uniprot_input
        st.session_state.pdb_data = None           # force reload
        st.session_state.protein_info = None
        st.session_state.mutation_result = None

    st.divider()
    st.header("Display Options")
    style = st.selectbox("Representation", ["cartoon", "stick", "sphere", "surface"])
    color = st.selectbox("Color Scheme",   ["confidence", "rainbow", "secondary"])

# ── Auto-fetch if we have an ID but no data yet ──────────────────────────────
if st.session_state.loaded_id and st.session_state.pdb_data is None:
    with st.spinner(f"Fetching {st.session_state.loaded_id} from AlphaFold..."):
        st.session_state.pdb_data    = fetch_pdb(st.session_state.loaded_id)
        st.session_state.protein_info = fetch_protein_info(st.session_state.loaded_id)

# ── Main content ─────────────────────────────────────────────────────────────
if st.session_state.pdb_data is None and not st.session_state.loaded_id:

    # Landing page
    st.info("👈 Select a protein from the sidebar or enter a UniProt ID to get started!")
    st.markdown("### 🧪 What is a UniProt ID?")
    st.markdown("""
    Every protein has a unique ID in the UniProt database — like a passport number.
    - `P69905` → Human Hemoglobin (carries oxygen in blood)
    - `P01308` → Human Insulin (regulates blood sugar)
    - `P38398` → BRCA1 (breast cancer related gene)
    """)

elif st.session_state.pdb_data is None:
    st.error(f"❌ Could not find protein '{st.session_state.loaded_id}'. Check the UniProt ID.")

else:
    pdb_data = st.session_state.pdb_data
    info     = st.session_state.protein_info

    # ── 3D Viewer + Info ─────────────────────────────────────────
    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("🧬 3D Structure")
        render_protein(pdb_data, style=style, color=color)
        if color == "confidence":
            st.markdown("""
            🔵 **90–100** Very high confidence &nbsp;|&nbsp;
            🟢 **70–90** Confident &nbsp;|&nbsp;
            🟡 **50–70** Low confidence &nbsp;|&nbsp;
            🔴 **<50** Very low confidence
            """)

    with col2:
        if info:
            st.subheader("📋 Protein Info")
            st.metric("Protein",   info["name"])
            st.metric("Organism",  info["organism"])
            st.metric("Length",    f"{info['length']} amino acids")
            st.markdown("**Function:**")
            st.info(info["function"][:400] + "..." if len(info["function"]) > 400 else info["function"])
            st.markdown(f"[UniProt Page](https://www.uniprot.org/uniprot/{st.session_state.loaded_id})")
            st.markdown(f"[AlphaFold Page](https://alphafold.ebi.ac.uk/entry/{st.session_state.loaded_id})")

    # ── pLDDT Chart ───────────────────────────────────────────────
    st.subheader("Confidence Score Per Residue (pLDDT)")
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
        bar_colors = [
            "#003f7f" if s >= 90 else
            "#1f9e89" if s >= 70 else
            "#f5a623" if s >= 50 else
            "#d73027" for s in scores
        ]
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=residues, y=scores,
            marker_color=bar_colors,
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

    # ── Mutation Lab ──────────────────────────────────────────────
    st.divider()
    st.subheader("🧪 Mutation Lab")
    st.markdown("Simulate swapping one amino acid and see how it affects stability")

    amino_acids = {
        "A - Alanine": "A",       "C - Cysteine": "C",
        "D - Aspartate": "D",     "E - Glutamate": "E",
        "F - Phenylalanine": "F", "G - Glycine": "G",
        "H - Histidine": "H",     "I - Isoleucine": "I",
        "K - Lysine": "K",        "L - Leucine": "L",
        "M - Methionine": "M",    "N - Asparagine": "N",
        "P - Proline": "P",       "Q - Glutamine": "Q",
        "R - Arginine": "R",      "S - Serine": "S",
        "T - Threonine": "T",     "V - Valine": "V",
        "W - Tryptophan": "W",    "Y - Tyrosine": "Y",
    }

    mut_col1, mut_col2, mut_col3 = st.columns([1, 1, 1])

    with mut_col1:
        position = st.number_input(
            "Position (residue number)",
            min_value=1,
            max_value=info["length"] if info else 9999,
            value=6,
            key="position_input"
        )
    with mut_col2:
        new_aa_label = st.selectbox(
            "Change to",
            list(amino_acids.keys()),
            key="new_aa_select"
        )
        new_aa = amino_acids[new_aa_label]
    with mut_col3:
        st.markdown("<br>", unsafe_allow_html=True)
        analyze_btn = st.button("🔬 Analyze Mutation", type="primary", use_container_width=True)

    if analyze_btn:
        original_aa = get_residue_at_position(pdb_data, position)

        if original_aa is None:
            st.error(f"❌ No residue found at position {position}.")
        elif original_aa == new_aa:
            st.warning(f"⚠️ Position {position} is already {new_aa}! Pick a different amino acid.")
        else:
            plddt       = get_plddt_at_position(pdb_data, position) or 80.0
            with st.spinner("⏳ Submitting to DynaMut2 server... (may take ~20 seconds)"):
                ddg_data = calculate_ddg(pdb_data, original_aa, new_aa, position, plddt)
            interpreted = interpret_ddg(ddg_data["ddg"], original_aa, new_aa, position)

            # Store result in session state so it persists
            st.session_state.mutation_result = {
                "position":    position,
                "original_aa": original_aa,
                "new_aa":      new_aa,
                "plddt":       plddt,
                "ddg_data":    ddg_data,
                "interpreted": interpreted,
            }

    # ── Show mutation result if it exists ─────────────────────────
    if st.session_state.mutation_result:
        r = st.session_state.mutation_result

        st.markdown(f"### Mutation: **{r['original_aa']}{r['position']}{r['new_aa']}**")

        res_col1, res_col2 = st.columns([2, 1])

        with res_col1:
            render_protein_with_mutation(pdb_data, r["position"])

        with res_col2:
            st.markdown("### Stability Analysis")
            st.metric(
                label="ΔΔG (kcal/mol)",
                value=f"{r['ddg_data']['ddg']:+.2f}",
                delta="destabilizing" if r['ddg_data']['ddg'] > 0.3 else "neutral",
                delta_color="inverse"
            )

            interp = r["interpreted"]
            st.markdown(
                f"<div style='background:{interp['color']}22;"
                f"border-left:4px solid {interp['color']};"
                f"padding:12px;border-radius:6px;margin:10px 0'>"
                f"<b>{interp['emoji']} {interp['level']}</b></div>",
                unsafe_allow_html=True
            )

            st.markdown("**What changed:**")
            st.markdown(f"- Hydrophobicity shift: `{r['ddg_data']['hydrophobicity_diff']}`")
            st.markdown(f"- Charge change: `{r['ddg_data']['charge_diff']}`")
            st.markdown(f"- Size change: `{r['ddg_data']['size_diff']}`")
            st.markdown(f"- Region confidence: `{int(r['plddt'])} pLDDT`")

        st.markdown("### Explanation")
        st.info(interp["explanation"])

        if r["original_aa"] == "E" and r["new_aa"] == "V" and r["position"] == 6:
            st.success(
                "🩸 **You just simulated the Sickle Cell Mutation!** "
                "This exact change (E6V in Hemoglobin) affects millions of people worldwide. "
                "One letter. One disease."
            )