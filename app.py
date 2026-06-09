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
from disease_lookup import get_disease_info
from ai_explainer import generate_disease_explanation
from ml_model.predict import predict_pathogenicity
from variant_tools import get_variant_scores
from protein_comparison import (
    extract_sequence_from_pdb,
    extract_plddt_from_pdb,
    align_sequences,
    compare_properties,
    classify_similarity,
    AA_PROPERTIES,
)
from conservation import get_conservation, interpret_conservation_for_mutation
from population_frequency import get_population_frequency

def _generate_pdf_report(history: list, df) -> bytes | None:
    """Generate a PDF report of all analyzed mutations."""
    try:
        from reportlab.lib.pagesizes import letter, A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.lib import colors
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer,
            Table, TableStyle, HRFlowable
        )
        from reportlab.lib.enums import TA_CENTER, TA_LEFT
        import io
        import datetime

        buffer = io.BytesIO()
        doc    = SimpleDocTemplate(
            buffer, pagesize=A4,
            topMargin=0.75*inch, bottomMargin=0.75*inch,
            leftMargin=0.75*inch, rightMargin=0.75*inch
        )

        styles = getSampleStyleSheet()
        story  = []

        # ── Title ─────────────────────────────────────────────
        title_style = ParagraphStyle(
            "Title",
            parent=styles["Title"],
            fontSize=20,
            textColor=colors.HexColor("#00d4ff"),
            spaceAfter=6,
        )
        story.append(Paragraph("🔬 Protein Mutation Analysis Report", title_style))
        story.append(Paragraph(
            f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
            styles["Normal"]
        ))
        story.append(HRFlowable(width="100%", thickness=1,
                                color=colors.HexColor("#00d4ff")))
        story.append(Spacer(1, 0.2*inch))

        # ── Summary ───────────────────────────────────────────
        story.append(Paragraph("Summary", styles["Heading2"]))
        summary_data = [
            ["Metric", "Value"],
            ["Total mutations analyzed", str(len(history))],
            ["Proteins studied",
             str(len(set(h["protein"] for h in history)))],
            ["Average ΔΔG",
             f"{sum(h['ddg'] for h in history)/len(history):+.2f} kcal/mol"],
            ["Likely harmful mutations",
             str(sum(1 for h in history
                     if "Destabilizing" in h.get("stability", "")))],
        ]
        summary_table = Table(summary_data, colWidths=[3*inch, 3*inch])
        summary_table.setStyle(TableStyle([
            ("BACKGROUND",  (0, 0), (-1, 0), colors.HexColor("#00d4ff")),
            ("TEXTCOLOR",   (0, 0), (-1, 0), colors.white),
            ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE",    (0, 0), (-1, -1), 10),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1),
             [colors.HexColor("#1a1a2e"), colors.HexColor("#0f0f1a")]),
            ("TEXTCOLOR",   (0, 1), (-1, -1), colors.white),
            ("GRID",        (0, 0), (-1, -1), 0.5, colors.HexColor("#333")),
            ("PADDING",     (0, 0), (-1, -1), 8),
        ]))
        story.append(summary_table)
        story.append(Spacer(1, 0.3*inch))

        # ── Mutations table ───────────────────────────────────
        story.append(Paragraph("Mutation Details", styles["Heading2"]))

        table_data = [["Mutation", "Protein", "ΔΔG", "Stability",
                       "Disease", "ClinVar", "ML Score"]]
        for h in history:
            table_data.append([
                h["mutation"],
                h["protein_name"][:20],
                f"{h['ddg']:+.2f}",
                h["stability"][:20],
                h.get("disease", "—")[:20],
                h.get("clinvar", "—")[:15],
                h.get("ml_score", "—"),
            ])

        mut_table = Table(
            table_data,
            colWidths=[0.9*inch, 1.5*inch, 0.7*inch, 1.5*inch,
                       1.5*inch, 1.0*inch, 0.7*inch]
        )
        mut_table.setStyle(TableStyle([
            ("BACKGROUND",  (0, 0), (-1, 0), colors.HexColor("#e74c3c")),
            ("TEXTCOLOR",   (0, 0), (-1, 0), colors.white),
            ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE",    (0, 0), (-1, -1), 8),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1),
             [colors.HexColor("#1a1a2e"), colors.HexColor("#0f0f1a")]),
            ("TEXTCOLOR",   (0, 1), (-1, -1), colors.white),
            ("GRID",        (0, 0), (-1, -1), 0.5, colors.HexColor("#333")),
            ("PADDING",     (0, 0), (-1, -1), 6),
            ("FONTNAME",    (0, 1), (-1, -1), "Helvetica"),
        ]))
        story.append(mut_table)
        story.append(Spacer(1, 0.3*inch))

        # ── Individual mutation details ───────────────────────
        story.append(Paragraph("Individual Mutation Analysis", styles["Heading2"]))

        for h in history:
            story.append(Paragraph(
                f"Mutation: {h['mutation']} in {h['protein_name']}",
                styles["Heading3"]
            ))
            details = [
                f"Position: {h['position']}  |  "
                f"Original: {h['original_aa']}  |  New: {h['new_aa']}",
                f"ΔΔG: {h['ddg']:+.2f} kcal/mol  |  "
                f"Stability: {h['stability']}",
                f"pLDDT at site: {h['plddt']}",
                f"Disease: {h.get('disease', 'Not analyzed')}",
                f"ClinVar: {h.get('clinvar', 'Not searched')}",
                f"ML Pathogenicity: {h.get('ml_score', 'Not analyzed')}",
            ]
            for detail in details:
                story.append(Paragraph(detail, styles["Normal"]))
            story.append(Spacer(1, 0.1*inch))
            story.append(HRFlowable(
                width="100%", thickness=0.5,
                color=colors.HexColor("#333")
            ))
            story.append(Spacer(1, 0.1*inch))

        # ── Footer ────────────────────────────────────────────
        story.append(Spacer(1, 0.3*inch))
        story.append(HRFlowable(
            width="100%", thickness=1,
            color=colors.HexColor("#00d4ff")
        ))
        story.append(Paragraph(
            "Generated by Protein Mutation Explorer · "
            "Data sources: AlphaFold, ClinVar, DynaMut2, AlphaMissense",
            ParagraphStyle("Footer", parent=styles["Normal"],
                           fontSize=8, textColor=colors.grey,
                           alignment=TA_CENTER)
        ))

        doc.build(story)
        return buffer.getvalue()

    except ImportError:
        st.error("Install reportlab: pip install reportlab")
        return None
    except Exception as e:
        st.error(f"PDF generation error: {e}")
        return None

# ── Page config ───────────────────────────────────────────────────
st.set_page_config(page_title="Protein Explorer", page_icon="🔬", layout="wide")

# ── Initialize session state ──────────────────────────────────────
if "pdb_data" not in st.session_state:
    st.session_state.pdb_data = None
if "protein_info" not in st.session_state:
    st.session_state.protein_info = None
if "loaded_id" not in st.session_state:
    st.session_state.loaded_id = ""
if "mutation_result" not in st.session_state:
    st.session_state.mutation_result = None
if "disease_result" not in st.session_state:
    st.session_state.disease_result = None
if "variant_scores" not in st.session_state:
    st.session_state.variant_scores = None
if "ml_result" not in st.session_state:
    st.session_state.ml_result = None
if "mutation_history" not in st.session_state:
    st.session_state.mutation_history = []
if "disease_history" not in st.session_state:
    st.session_state.disease_history = {}
if "pdb_data2" not in st.session_state:
    st.session_state.pdb_data2 = None
if "protein_info2" not in st.session_state:
    st.session_state.protein_info2 = None
if "loaded_id2" not in st.session_state:
    st.session_state.loaded_id2 = ""
if "comparison_result" not in st.session_state:
    st.session_state.comparison_result = None
if "conservation" not in st.session_state:
    st.session_state.conservation = None
if "population_freq" not in st.session_state:
    st.session_state.population_freq = None
# ── Header ────────────────────────────────────────────────────────
st.title("🔬 Protein Structure Explorer")

# ── Sidebar ───────────────────────────────────────────────────────
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
            st.session_state.loaded_id    = uid
            st.session_state.pdb_data     = None
            st.session_state.protein_info = None
            st.session_state.mutation_result = None
            st.session_state.disease_result  = None
            st.session_state.variant_scores  = None

    st.divider()

    uniprot_input = st.text_input(
        "Or enter UniProt ID manually:",
        value=st.session_state.loaded_id
    ).strip().upper()

    if st.button("🔍 Load Protein", type="primary", use_container_width=True):
        st.session_state.loaded_id    = uniprot_input
        st.session_state.pdb_data     = None
        st.session_state.protein_info = None
        st.session_state.mutation_result = None
        st.session_state.disease_result  = None
        st.session_state.variant_scores  = None

    st.divider()
    st.header("🎨 Display Options")
    style = st.selectbox("Representation", ["cartoon", "stick", "sphere", "surface"])
    color = st.selectbox("Color Scheme",   ["confidence", "rainbow", "secondary"])

# ── PASTE THE ENTIRE MAIN CONTENT BLOCK HERE ─────────────────────
# (everything starting from "# ── Auto-fetch" that I gave you)

# ── Auto-fetch ────────────────────────────────────────────────────
if st.session_state.loaded_id and st.session_state.pdb_data is None:
    with st.spinner(f"Fetching {st.session_state.loaded_id} from AlphaFold..."):
        st.session_state.pdb_data     = fetch_pdb(st.session_state.loaded_id)
        st.session_state.protein_info = fetch_protein_info(st.session_state.loaded_id)

# ── Landing page ──────────────────────────────────────────────────
if st.session_state.pdb_data is None and not st.session_state.loaded_id:
    st.markdown("""
    <div style='text-align:center; padding:60px 20px'>
        <h1 style='font-size:3em'>🔬</h1>
        <h2>Protein Mutation Explorer</h2>
        <p style='color:#aaa; font-size:1.1em'>
            Visualize proteins in 3D · Simulate mutations · Predict diseases with AI
        </p>
    </div>
    """, unsafe_allow_html=True)

    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("""
        <div style='background:#1a1a2e; border-left:4px solid #00d4ff;
        padding:20px; border-radius:8px; text-align:center'>
        <h3>🧬 3D Structure</h3>
        <p style='color:#aaa'>Fetch any protein from DeepMind's AlphaFold
        database and explore it interactively in 3D</p>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown("""
        <div style='background:#1a1a2e; border-left:4px solid #ff6b6b;
        padding:20px; border-radius:8px; text-align:center'>
        <h3>🧪 Mutation Lab</h3>
        <p style='color:#aaa'>Simulate amino acid changes and calculate
        structural stability impact using DynaMut2</p>
        </div>""", unsafe_allow_html=True)
    with c3:
        st.markdown("""
        <div style='background:#1a1a2e; border-left:4px solid #2ecc71;
        padding:20px; border-radius:8px; text-align:center'>
        <h3>🏥 Disease Prediction</h3>
        <p style='color:#aaa'>Predict disease impact using XGBoost,
        ClinVar, SIFT, PolyPhen, AlphaMissense + AI explanation</p>
        </div>""", unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### 🚀 Quick Start — Try These")
    q1, q2, q3, q4 = st.columns(4)
    with q1:
        if st.button("🩸 Sickle Cell\nHemoglobin", use_container_width=True):
            st.session_state.loaded_id    = "P69905"
            st.session_state.pdb_data     = None
            st.session_state.protein_info = None
            st.rerun()
    with q2:
        if st.button("🎗️ BRCA1\nBreast Cancer", use_container_width=True):
            st.session_state.loaded_id    = "P38398"
            st.session_state.pdb_data     = None
            st.session_state.protein_info = None
            st.rerun()
    with q3:
        if st.button("💉 Insulin\nDiabetes", use_container_width=True):
            st.session_state.loaded_id    = "P01308"
            st.session_state.pdb_data     = None
            st.session_state.protein_info = None
            st.rerun()
    with q4:
        if st.button("🦠 COVID-19\nSpike Protein", use_container_width=True):
            st.session_state.loaded_id    = "P0DTC2"
            st.session_state.pdb_data     = None
            st.session_state.protein_info = None
            st.rerun()

elif st.session_state.pdb_data is None:
    st.error(f"❌ Could not find protein '{st.session_state.loaded_id}'. Check the UniProt ID.")

else:
    pdb_data = st.session_state.pdb_data
    info     = st.session_state.protein_info

    # ── Protein header ────────────────────────────────────────────
    protein_name = info["name"] if info else st.session_state.loaded_id
    gene_symbol  = info.get("gene_symbol", "") if info else ""
    organism     = info.get("organism", "") if info else ""

    st.markdown(
        f"<div style='background:#1a1a2e; padding:16px 20px; "
        f"border-radius:8px; margin-bottom:16px; "
        f"border-left:4px solid #00d4ff'>"
        f"<h3 style='margin:0; color:#00d4ff'>{protein_name}</h3>"
        f"<p style='margin:4px 0 0 0; color:#aaa'>"
        f"Gene: <b style='color:white'>{gene_symbol}</b> &nbsp;|&nbsp; "
        f"Organism: <b style='color:white'>{organism}</b> &nbsp;|&nbsp; "
        f"UniProt: <b style='color:white'>{st.session_state.loaded_id}</b>"
        f"</p></div>",
        unsafe_allow_html=True
    )

    # ── TABS ──────────────────────────────────────────────────────
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "🧬 Structure",
        "🧪 Mutation Lab",
        "🏥 Disease Prediction",
        "📊 History",
        "🔬 Compare Proteins",
        "ℹ️ How to Use",
    ])

    # ════════════════════════════════════════════════════════════
    # TAB 1 — Structure
    # ════════════════════════════════════════════════════════════
    with tab1:
        col1, col2 = st.columns([2, 1])

        with col1:
            st.subheader("🧬 3D Structure")
            render_protein(pdb_data, style=style, color=color)
            if color == "confidence":
                st.markdown("""
                🔵 **90–100** Very high &nbsp;|&nbsp;
                🟢 **70–90** Confident &nbsp;|&nbsp;
                🟡 **50–70** Low &nbsp;|&nbsp;
                🔴 **<50** Very low confidence
                """)

        with col2:
            if info:
                st.subheader("📋 Protein Info")
                st.metric("Length", f"{info['length']} amino acids")
                st.markdown("**Function:**")
                fn = info["function"]
                st.info(fn[:400] + "..." if len(fn) > 400 else fn)
                st.markdown(f"[UniProt ↗](https://www.uniprot.org/uniprot/{st.session_state.loaded_id})")
                st.markdown(f"[AlphaFold ↗](https://alphafold.ebi.ac.uk/entry/{st.session_state.loaded_id})")

        # pLDDT chart
        st.subheader("📊 Confidence Score Per Residue (pLDDT)")
        residues, pscores = [], []
        seen = set()
        for line in pdb_data.split("\n"):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_num = int(line[22:26].strip())
                if res_num not in seen:
                    try:
                        residues.append(res_num)
                        pscores.append(float(line[60:66].strip()))
                        seen.add(res_num)
                    except:
                        pass

        if residues:
            bar_colors = [
                "#003f7f" if s >= 90 else
                "#1f9e89" if s >= 70 else
                "#f5a623" if s >= 50 else
                "#d73027" for s in pscores
            ]
            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=residues, y=pscores,
                marker_color=bar_colors,
                hovertemplate="Residue %{x}<br>pLDDT: %{y:.1f}<extra></extra>"
            ))
            fig.add_hline(y=90, line_dash="dot", line_color="#003f7f",
                          annotation_text="Very High (90)")
            fig.add_hline(y=70, line_dash="dot", line_color="#1f9e89",
                          annotation_text="Confident (70)")
            fig.add_hline(y=50, line_dash="dot", line_color="#f5a623",
                          annotation_text="Low (50)")
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

    # ════════════════════════════════════════════════════════════
    # TAB 2 — Mutation Lab
    # ════════════════════════════════════════════════════════════
    with tab2:
        st.subheader("🧪 Mutation Lab")
        st.markdown("Simulate swapping one amino acid and see how it affects protein stability")

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
                value=6, key="position_input"
            )
        with mut_col2:
            new_aa_label = st.selectbox(
                "Change to", list(amino_acids.keys()),
                key="new_aa_select"
            )
            new_aa = amino_acids[new_aa_label]
        with mut_col3:
            st.markdown("<br>", unsafe_allow_html=True)
            analyze_btn = st.button(
                "🔬 Analyze Mutation",
                type="primary",
                use_container_width=True
            )

        if analyze_btn:
            original_aa = get_residue_at_position(pdb_data, position)
            if original_aa is None:
                st.error(f"❌ No residue found at position {position}.")
            elif original_aa == new_aa:
                st.warning(f"⚠️ Position {position} is already {new_aa}!")
            else:
                plddt = get_plddt_at_position(pdb_data, position) or 80.0
                with st.spinner("⏳ Calculating stability... (may take ~20 seconds)"):
                    ddg_data = calculate_ddg(pdb_data, original_aa, new_aa, position, plddt)
                interpreted = interpret_ddg(ddg_data["ddg"], original_aa, new_aa, position)
                st.session_state.mutation_result = {
                    "position":    position,
                    "original_aa": original_aa,
                    "new_aa":      new_aa,
                    "plddt":       plddt,
                    "ddg_data":    ddg_data,
                    "interpreted": interpreted,
                }
                st.session_state.disease_result  = None
                st.session_state.variant_scores  = None

                # ── Save to history ───────────────────────────────────
                import datetime

                history_entry = {
                    "mutation": f"{original_aa}{position}{new_aa}",
                    "position": position,
                    "original_aa": original_aa,
                    "new_aa": new_aa,
                    "ddg": round(ddg_data["ddg"], 2),
                    "stability": interpreted["level"],
                    "plddt": round(plddt, 1),
                    "protein": st.session_state.loaded_id,
                    "protein_name": info["name"][:30] if info else "Unknown",
                    "timestamp": datetime.datetime.now().strftime("%H:%M:%S"),
                }
                existing = [h["mutation"] for h in st.session_state.mutation_history]
                if history_entry["mutation"] not in existing:
                    st.session_state.mutation_history.append(history_entry)

        if st.session_state.mutation_result:
            r = st.session_state.mutation_result
            st.markdown(f"### Result: **{r['original_aa']}{r['position']}{r['new_aa']}**")

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
                st.markdown(f"- Hydrophobicity: `{r['ddg_data']['hydrophobicity_diff']}`")
                st.markdown(f"- Charge: `{r['ddg_data']['charge_diff']}`")
                st.markdown(f"- Size: `{r['ddg_data']['size_diff']}`")
                st.markdown(f"- pLDDT: `{int(r['plddt'])}/100`")

            st.markdown("### 📝 Explanation")
            st.info(interp["explanation"])

            if r["original_aa"] == "E" and r["new_aa"] == "V" and r["position"] == 6:
                st.success(
                    "🩸 **You just simulated the Sickle Cell Mutation (E6V)!** "
                    "One amino acid change. Millions of lives affected."
                )

            st.info("👉 Go to the **Disease Prediction** tab to see what disease this causes!")

        else:
            st.markdown("""
            <div style='text-align:center; padding:40px; color:#aaa'>
                <h3>👆 Enter a position and amino acid above</h3>
                <p>Try position <b>6</b> → <b>V - Valine</b> on Hemoglobin
                to simulate the Sickle Cell mutation</p>
            </div>
            """, unsafe_allow_html=True)

    # ════════════════════════════════════════════════════════════
    # TAB 3 — Disease Prediction
    # ════════════════════════════════════════════════════════════
    with tab3:
        st.subheader("🏥 Disease Prediction")

        if not st.session_state.mutation_result:
            st.markdown("""
            <div style='text-align:center; padding:40px; color:#aaa;
            background:#1a1a2e; border-radius:8px'>
                <h3>⚠️ No mutation analyzed yet</h3>
                <p>Go to the <b>Mutation Lab</b> tab first,
                analyze a mutation, then come back here.</p>
            </div>
            """, unsafe_allow_html=True)
        else:
            r = st.session_state.mutation_result
            st.markdown(
                f"Predicting disease impact of mutation "
                f"**{r['original_aa']}{r['position']}{r['new_aa']}**"
            )

            dp_col1, dp_col2 = st.columns([1, 1])
            with dp_col1:
                auto_gene  = info.get("gene_symbol", "Unknown") if info else "Unknown"
                gene_input = st.text_input(
                    "Gene symbol",
                    value=auto_gene,
                    help="e.g. HBB, BRCA1, TP53, CFTR",
                    key="gene_input"
                )
            with dp_col2:
                st.markdown("<br>", unsafe_allow_html=True)
                predict_btn = st.button(
                    "🔬 Predict Disease",
                    type="primary",
                    use_container_width=True
                )

            if predict_btn:
                with st.spinner("🧬 Running ML model..."):
                    ml_result = predict_pathogenicity(
                        r["original_aa"], r["new_aa"],
                        r["plddt"], gene=gene_input
                    )
                st.session_state.ml_result = ml_result

                with st.spinner("🔍 Searching ClinVar..."):
                    disease_info = get_disease_info(
                        gene_input, r["original_aa"],
                        r["position"], r["new_aa"]
                    )

                with st.spinner("🔬 Running SIFT, PolyPhen & AlphaMissense..."):
                    variant_scores = get_variant_scores(
                        gene=gene_input,
                        uniprot_id=st.session_state.loaded_id,
                        original_aa=r["original_aa"],
                        position=r["position"],
                        new_aa=r["new_aa"]
                    )
                st.session_state.variant_scores = variant_scores

                with st.spinner("🤖 Generating AI explanation..."):
                    ai_result = generate_disease_explanation(
                        gene=gene_input,
                        mutation=f"{r['original_aa']}{r['position']}{r['new_aa']}",
                        clinvar_data=disease_info.get("clinvar"),
                        omim_data=disease_info.get("omim"),
                        ddg_data=r["ddg_data"],
                        ml_score=ml_result["score"],
                        plddt=r["plddt"],
                    )

                # Save disease result
                st.session_state.disease_result = {
                    "ml": ml_result,
                    "disease": disease_info,
                    "ai": ai_result,
                }

                # ← Attach disease info to history entry
                mutation_key = f"{r['original_aa']}{r['position']}{r['new_aa']}"
                for entry in st.session_state.mutation_history:
                    if entry["mutation"] == mutation_key:
                        entry["disease"] = ai_result.get("disease_name", "Unknown")
                        entry["risk"] = ai_result.get("risk_level", "Unknown")
                        entry["clinvar"] = disease_info.get("clinvar", {}).get("significance", "Not found")
                        entry["ml_score"] = f"{ml_result['score']:.0%}"
                        entry["gene"] = gene_input
                        break

            if st.session_state.disease_result:
                dr = st.session_state.disease_result
                ml = dr["ml"]
                ai = dr["ai"]
                cv = dr["disease"].get("clinvar", {})
                om = dr["disease"].get("omim", {})

                # Top summary banner
                st.markdown(
                    f"<div style='background:{ai['risk_color']}22;"
                    f"border:2px solid {ai['risk_color']};"
                    f"border-radius:12px;padding:20px;margin:10px 0'>"
                    f"<h2 style='margin:0;color:{ai['risk_color']}'>"
                    f"{ai['risk_emoji']} {ai['disease_name']}</h2>"
                    f"<p style='margin:5px 0 0 0;color:#aaa'>"
                    f"Risk: <b style='color:{ai['risk_color']}'>{ai['risk_level']}</b>"
                    f" &nbsp;|&nbsp; Severity: <b>{ai['severity']}</b>"
                    f" &nbsp;|&nbsp; System: <b>{ai['body_system']}</b></p>"
                    f"</div>",
                    unsafe_allow_html=True
                )

                # Three source columns
                s1, s2, s3 = st.columns(3)
                with s1:
                    st.markdown("#### 🤖 ML Model")
                    xgb_score  = ml.get("xgb_score", ml["score"])
                    gene_known = ml.get("gene_known", False)
                    st.markdown(
                        f"<div style='background:{ml['color']}22;"
                        f"border-left:4px solid {ml['color']};"
                        f"padding:12px;border-radius:6px'>"
                        f"<b style='font-size:1.4em'>{ml['score']:.0%}</b><br>"
                        f"{ml['emoji']} {ml['classification']}<br>"
                        f"<small style='color:#aaa'>{ml['confidence']}</small><br><br>"
                        f"<small>XGBoost: {xgb_score:.0%}</small><br>"
                        f"<small>Gene: {'✅ Known' if gene_known else '❓ Unknown'}</small>"
                        f"</div>",
                        unsafe_allow_html=True
                    )


                with s2:
                    st.markdown("#### 📋 ClinVar")
                    if cv and cv.get("status") == "found":
                        sc = cv.get("color", "#95a5a6")
                        st.markdown(
                            f"<div style='background:{sc}22;"
                            f"border-left:4px solid {sc};"
                            f"padding:12px;border-radius:6px'>"
                            f"{cv['emoji']} <b>{cv['significance']}</b><br>"
                            f"<small>{cv['review_status']}</small><br>"
                            f"<small style='color:#aaa'>"
                            f"{', '.join(cv['diseases'][:2])}</small>"
                            f"</div>",
                            unsafe_allow_html=True
                        )
                        if cv.get("clinvar_url"):
                            st.markdown(f"[View on ClinVar ↗]({cv['clinvar_url']})")
                    else:
                        st.markdown(
                            "<div style='background:#95a5a622;"
                            "border-left:4px solid #95a5a6;"
                            "padding:12px;border-radius:6px'>"
                            "❓ Not in ClinVar<br>"
                            "<small>Novel or unstudied mutation</small>"
                            "</div>", unsafe_allow_html=True
                        )

                with s3:
                    st.markdown("#### 🧬 Inheritance")
                    inheritance = (
                        ml.get("inheritance") or
                        ai.get("inheritance") or
                        cv.get("inheritance", "Not determined")
                    )
                    if inheritance == "Unknown":
                        inheritance = "Not determined"
                    st.markdown(
                        f"<div style='background:#3498db22;"
                        f"border-left:4px solid #3498db;"
                        f"padding:12px;border-radius:6px'>"
                        f"<b>{inheritance}</b><br><br>"
                        f"<small style='color:#aaa'>Key symptoms:</small><br>"
                        + "".join([f"<small>• {s}</small><br>"
                                   for s in ai['key_symptoms'][:4]]) +
                        f"</div>", unsafe_allow_html=True
                    )

                # Variant scoring tools
                if st.session_state.variant_scores:
                    vs = st.session_state.variant_scores
                    sp_vs = vs.get("sift_polyphen", {})
                    am_vs = vs.get("alphamissense")
                    ov_vs = vs.get("overall_verdict", {})

                    st.markdown("---")
                    st.markdown("#### 🔬 Variant Scoring Tools")
                    st.caption("Industry-standard tools used in clinical genetics labs worldwide")

                    t1, t2, t3, t4 = st.columns(4)
                    with t1:
                        st.markdown("**SIFT**")
                        sl = sp_vs.get("sift_prediction","Not available").replace("_"," ").title()
                        badge = " · 📚" if sp_vs.get("source") == "curated" else " · 🌐"
                        st.markdown(
                            f"<div style='background:{sp_vs.get('sift_color','#95a5a6')}22;"
                            f"border-left:4px solid {sp_vs.get('sift_color','#95a5a6')};"
                            f"padding:10px;border-radius:6px'>"
                            f"{sp_vs.get('sift_emoji','⚪')} <b>{sl}</b>"
                            + (f"<br><small>Score: {sp_vs['sift_score']}</small>"
                               if sp_vs.get("sift_score") is not None else "") +
                            f"<br><small style='color:#aaa'>Evolutionary{badge}</small>"
                            f"</div>", unsafe_allow_html=True
                        )
                    with t2:
                        st.markdown("**PolyPhen-2**")
                        pl = sp_vs.get("polyphen_prediction","Not available").replace("_"," ").title()
                        st.markdown(
                            f"<div style='background:{sp_vs.get('polyphen_color','#95a5a6')}22;"
                            f"border-left:4px solid {sp_vs.get('polyphen_color','#95a5a6')};"
                            f"padding:10px;border-radius:6px'>"
                            f"{sp_vs.get('polyphen_emoji','⚪')} <b>{pl}</b>"
                            + (f"<br><small>Score: {sp_vs['polyphen_score']}</small>"
                               if sp_vs.get("polyphen_score") is not None else "") +
                            f"<br><small style='color:#aaa'>Structure + sequence</small>"
                            f"</div>", unsafe_allow_html=True
                        )
                    with t3:
                        st.markdown("**AlphaMissense**")
                        if am_vs and am_vs.get("status") == "found":
                            st.markdown(
                                f"<div style='background:{am_vs['color']}22;"
                                f"border-left:4px solid {am_vs['color']};"
                                f"padding:10px;border-radius:6px'>"
                                f"{am_vs['emoji']} <b>{am_vs['classification']}</b>"
                                f"<br><small>Score: {am_vs['score']}</small>"
                                f"<br><small style='color:#aaa'>DeepMind AI</small>"
                                f"</div>", unsafe_allow_html=True
                            )
                        else:
                            st.markdown(
                                "<div style='background:#95a5a622;"
                                "border-left:4px solid #95a5a6;"
                                "padding:10px;border-radius:6px'>"
                                "⚪ <b>Not available</b>"
                                "<br><small style='color:#aaa'>DeepMind AI</small>"
                                "</div>", unsafe_allow_html=True
                            )
                    with t4:
                        st.markdown("**Overall Verdict**")
                        st.markdown(
                            f"<div style='background:{ov_vs.get('color','#95a5a6')}22;"
                            f"border:2px solid {ov_vs.get('color','#95a5a6')};"
                            f"padding:10px;border-radius:6px;text-align:center'>"
                            f"<b style='font-size:1.2em'>"
                            f"{ov_vs.get('emoji','⚪')} {ov_vs.get('verdict','Unknown')}</b>"
                            f"<br><small style='color:#aaa'>Combined</small>"
                            f"</div>", unsafe_allow_html=True
                        )

                    with st.expander("ℹ️ What do these tools mean?"):
                        st.markdown("""
                        **SIFT** — evolutionary conservation. If no species ever had this
                        change → dangerous. Least reliable for structural proteins.

                        **PolyPhen-2** — 3D structure + sequence. More reliable than SIFT
                        for repeat-domain proteins like collagen.

                        **AlphaMissense** — DeepMind's 2023 AI, trained on 71M mutations.
                        Highest weight in the overall verdict.

                        **Weighting:** AlphaMissense ×3 · PolyPhen-2 ×2 · SIFT ×1
                        """)

                # Combined verdict
                st.markdown("---")
                st.markdown("#### 🎯 Combined Verdict")
                signals = []
                if cv and cv.get("status") == "found":
                    sig = cv.get("level", "unknown")
                    if sig in ["pathogenic", "likely_pathogenic"]:
                        signals.append({"source": "ClinVar", "score": 1.0, "weight": 3})
                    elif sig in ["benign", "likely_benign"]:
                        signals.append({"source": "ClinVar", "score": 0.0, "weight": 3})
                    else:
                        signals.append({"source": "ClinVar", "score": 0.5, "weight": 1})
                signals.append({"source": "ML Model", "score": ml["score"], "weight": 2})
                ddg = r["ddg_data"].get("ddg", 0)
                ddg_score = 0.9 if ddg > 1.5 else 0.7 if ddg > 0.8 else 0.5 if ddg > 0.3 else 0.2
                signals.append({"source": "ΔΔG Stability", "score": ddg_score, "weight": 1})

                tw = sum(s["weight"] for s in signals)
                ws = sum(s["score"] * s["weight"] for s in signals)
                fs = ws / tw if tw > 0 else 0.5

                if fs >= 0.75:
                    vd, vc, ve = "Pathogenic",           "#e74c3c", "🔴"
                    vt = "Multiple sources agree this mutation is disease-causing."
                elif fs >= 0.55:
                    vd, vc, ve = "Likely Pathogenic",    "#e67e22", "🟠"
                    vt = "Evidence suggests this mutation is likely harmful."
                elif fs >= 0.40:
                    vd, vc, ve = "Uncertain Significance","#95a5a6", "⚪"
                    vt = "Evidence is mixed. Further analysis recommended."
                elif fs >= 0.25:
                    vd, vc, ve = "Likely Benign",        "#f1c40f", "🟡"
                    vt = "Evidence suggests this mutation is likely harmless."
                else:
                    vd, vc, ve = "Benign",               "#2ecc71", "🟢"
                    vt = "Multiple sources agree this mutation is harmless."

                st.markdown(
                    f"<div style='background:{vc}22;border:2px solid {vc};"
                    f"border-radius:12px;padding:20px;margin:10px 0'>"
                    f"<h3 style='margin:0;color:{vc}'>{ve} {vd}</h3>"
                    f"<p style='margin:8px 0 0 0;color:#ccc'>{vt}</p>"
                    f"</div>", unsafe_allow_html=True
                )

                v1, v2, v3 = st.columns(3)
                for i, sig in enumerate(signals[:3]):
                    pct   = sig["score"] * 100
                    col   = "#e74c3c" if pct >= 75 else "#e67e22" if pct >= 55 else "#95a5a6" if pct >= 40 else "#2ecc71"
                    with [v1, v2, v3][i]:
                        st.markdown(
                            f"<div style='background:#1a1a2e;border-left:3px solid {col};"
                            f"padding:8px;border-radius:4px;text-align:center'>"
                            f"<small style='color:#aaa'>{sig['source']}</small><br>"
                            f"<b style='color:{col}'>{pct:.0f}%</b>"
                            f"</div>", unsafe_allow_html=True
                        )

                # AI explanation
                st.markdown("---")
                st.markdown("#### 🗣️ AI Explanation")
                st.info(ai["explanation"])

                if om and om.get("description"):
                    with st.expander("📚 Full Disease Description"):
                        st.markdown(om["description"])
                        if om.get("omim_url"):
                            st.markdown(f"[View on OMIM ↗]({om['omim_url']})")

                if (r["original_aa"] == "E" and r["new_aa"] == "V"
                        and r["position"] == 6):
                    st.success(
                        "🩸 **You just analysed the Sickle Cell Mutation (E6V)!** "
                        "One amino acid change affects ~100 million people worldwide."
                    )
    # ════════════════════════════════════════════════════════════
    # TAB 4 — History & Export
    # ════════════════════════════════════════════════════════════
    with tab4:
        st.subheader("📊 Mutation History & Export")

        if not st.session_state.mutation_history:
            st.markdown("""
            <div style='text-align:center; padding:40px; color:#aaa;
            background:#1a1a2e; border-radius:8px'>
                <h3>📭 No mutations analyzed yet</h3>
                <p>Go to <b>Mutation Lab</b>, analyze some mutations,
                then come back here to compare them.</p>
            </div>
            """, unsafe_allow_html=True)

        else:
            history = st.session_state.mutation_history
            # ── Summary metrics ───────────────────────────────────
            m1, m2, m3, m4 = st.columns(4)
            with m1:
                st.metric("Mutations Analyzed", len(history))
            with m2:
                pathogenic_count = sum(
                    1 for h in history
                    if "Destabilizing" in h.get("stability", "")
                    or h.get("risk") in ["High", "Very High"]
                )
                st.metric("Likely Harmful", pathogenic_count)
            with m3:
                avg_ddg = sum(h["ddg"] for h in history) / len(history)
                st.metric("Avg ΔΔG", f"{avg_ddg:+.2f} kcal/mol")
            with m4:
                proteins = len(set(h["protein"] for h in history))
                st.metric("Proteins Studied", proteins)

            st.markdown("---")

            # ── Comparison table ──────────────────────────────────
            st.markdown("#### 🔬 All Mutations Compared")

            import pandas as pd

            # Build display dataframe
            rows = []
            for h in history:
                # Color code the stability
                stability_emoji = (
                    "🔴" if "Highly" in h.get("stability", "") else
                    "🟠" if "Moderately" in h.get("stability", "") else
                    "🟡" if "Mildly" in h.get("stability", "") else
                    "🟢"
                )
                rows.append({
                    "Mutation": h["mutation"],
                    "Protein": h["protein_name"],
                    "Position": h["position"],
                    "ΔΔG": f"{h['ddg']:+.2f}",
                    "Stability": f"{stability_emoji} {h['stability']}",
                    "pLDDT": h["plddt"],
                    "Disease": h.get("disease", "—"),
                    "ClinVar": h.get("clinvar", "—"),
                    "ML Score": h.get("ml_score", "—"),
                    "Risk": h.get("risk", "—"),
                    "Time": h["timestamp"],
                })
            df = pd.DataFrame(rows)
            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "ΔΔG": st.column_config.TextColumn("ΔΔG (kcal/mol)"),
                    "pLDDT": st.column_config.ProgressColumn(
                        "pLDDT", min_value=0, max_value=100, format="%d"
                    ),
                }
            )
            # ── ΔΔG comparison chart ──────────────────────────────
            st.markdown("#### 📈 ΔΔG Stability Comparison")
            st.caption("Higher bars = more destabilizing = more likely harmful")

            mutations = [h["mutation"] for h in history]
            ddg_vals = [h["ddg"] for h in history]
            bar_cols = [
                "#e74c3c" if d > 1.5 else
                "#e67e22" if d > 0.8 else
                "#f1c40f" if d > 0.3 else
                "#2ecc71" for d in ddg_vals
            ]
            fig = go.Figure()
            fig.add_trace(go.Bar(
                x=mutations,
                y=ddg_vals,
                marker_color=bar_cols,
                text=[f"{d:+.2f}" for d in ddg_vals],
                textposition="outside",
                hovertemplate="<b>%{x}</b><br>ΔΔG: %{y:+.2f} kcal/mol<extra></extra>"
            ))
            fig.add_hline(
                y=1.5, line_dash="dot", line_color="#e74c3c",
                annotation_text="Highly destabilizing threshold"
            )
            fig.add_hline(
                y=0.3, line_dash="dot", line_color="#2ecc71",
                annotation_text="Neutral threshold"
            )
            fig.update_layout(
                height=350,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis_title="Mutation",
                yaxis_title="ΔΔG (kcal/mol)",
                showlegend=False,
                font=dict(color="white"),
            )
            st.plotly_chart(fig, use_container_width=True)
            # ── pLDDT chart ───────────────────────────────────────
            st.markdown("#### 🎯 Mutation Site Confidence (pLDDT)")
            st.caption("Higher = mutation site is in a well-structured region")

            fig2 = go.Figure()
            fig2.add_trace(go.Bar(
                x=mutations,
                y=[h["plddt"] for h in history],
                marker_color=[
                    "#003f7f" if h["plddt"] >= 90 else
                    "#1f9e89" if h["plddt"] >= 70 else
                    "#f5a623" if h["plddt"] >= 50 else
                    "#d73027" for h in history
                ],
                text=[f"{h['plddt']:.0f}" for h in history],
                textposition="outside",
                hovertemplate="<b>%{x}</b><br>pLDDT: %{y:.1f}<extra></extra>"
            ))
            fig2.update_layout(
                height=300,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis_title="Mutation",
                yaxis_title="pLDDT Score",
                yaxis_range=[0, 110],
                showlegend=False,
                font=dict(color="white"),
            )
            st.plotly_chart(fig2, use_container_width=True)
            # ── Export section ────────────────────────────────────
            st.markdown("---")
            st.markdown("#### 💾 Export Results")

            ex1, ex2, ex3 = st.columns(3)

            # CSV Export
            with ex1:
                csv_data = df.to_csv(index=False)
                st.download_button(
                    label="📥 Download CSV",
                    data=csv_data,
                    file_name="mutation_analysis.csv",
                    mime="text/csv",
                    use_container_width=True,
                )
            # JSON Export
            with ex2:
                import json

                json_data = json.dumps(history, indent=2)
                st.download_button(
                    label="📥 Download JSON",
                    data=json_data,
                    file_name="mutation_analysis.json",
                    mime="application/json",
                    use_container_width=True,
                )
            # PDF Report Export
            with ex3:
                if st.button("📄 Generate PDF Report",
                             use_container_width=True,
                             type="primary"):
                    pdf_bytes = _generate_pdf_report(
                        st.session_state.mutation_history, df
                    )
                    if pdf_bytes:
                        st.download_button(
                            label="📥 Download PDF",
                            data=pdf_bytes,
                            file_name="mutation_report.pdf",
                            mime="application/pdf",
                            use_container_width=True,
                        )

            # Clear history
            st.markdown("---")
            if st.button("🗑️ Clear History", type="secondary"):
                st.session_state.mutation_history = []
                st.session_state.disease_history = {}
                st.rerun()

    # ════════════════════════════════════════════════════════════
    # TAB 5 — Compare Proteins
    # ════════════════════════════════════════════════════════════
    with tab5:
        st.subheader("🔬 Protein Comparison")
        st.markdown("Compare two proteins side by side — structure, sequence, and properties")

        # ── Protein selectors ─────────────────────────────────
        st.markdown("#### Select Two Proteins to Compare")

        cp1, cp2 = st.columns(2)

        with cp1:
            st.markdown(
                f"<div style='background:#1a1a2e; border-left:4px solid #00d4ff;"
                f"padding:12px; border-radius:8px'>"
                f"<b style='color:#00d4ff'>Protein A</b><br>"
                f"<small style='color:#aaa'>Currently loaded: "
                f"{st.session_state.loaded_id}</small></div>",
                unsafe_allow_html=True
            )
            st.caption("Protein A is your currently loaded protein")

        with cp2:
            st.markdown(
                "<div style='background:#1a1a2e; border-left:4px solid #ff6b6b;"
                "padding:12px; border-radius:8px'>"
                "<b style='color:#ff6b6b'>Protein B</b></div>",
                unsafe_allow_html=True
            )
            comp_examples = {
                "Select...":                "",
                "HBB — Beta Hemoglobin":    "P68871",
                "HBA1 — Alpha Hemoglobin":  "P69905",
                "Myoglobin":                "P02144",
                "Insulin":                  "P01308",
                "BRCA2":                    "P51587",
                "TP53":                     "P04637",
                "KRAS":                     "P01116",
                "Collagen II":              "P02458",
                "ACE2 (COVID receptor)":    "Q9BYF1",
            }

            comp_select = st.selectbox(
                "Quick select Protein B:",
                list(comp_examples.keys()),
                key="comp_select"
            )

            comp_manual = st.text_input(
                "Or enter UniProt ID:",
                value=comp_examples.get(comp_select, ""),
                key="comp_manual"
            ).strip().upper()

            compare_btn = st.button(
                "🔬 Compare Proteins",
                type="primary",
                use_container_width=True
            )

        # ── Run comparison ────────────────────────────────────
        if compare_btn and comp_manual:
            if comp_manual == st.session_state.loaded_id:
                st.warning("⚠️ Select a different protein to compare!")
            else:
                with st.spinner(f"Fetching Protein B ({comp_manual})..."):
                    pdb2  = fetch_pdb(comp_manual)
                    info2 = fetch_protein_info(comp_manual)

                if not pdb2:
                    st.error(f"❌ Could not find protein '{comp_manual}'")
                else:
                    st.session_state.pdb_data2     = pdb2
                    st.session_state.protein_info2 = info2
                    st.session_state.loaded_id2    = comp_manual

                    with st.spinner("🧬 Aligning sequences..."):
                        seq1   = extract_sequence_from_pdb(pdb_data)
                        seq2   = extract_sequence_from_pdb(pdb2)
                        plddt1 = extract_plddt_from_pdb(pdb_data)
                        plddt2 = extract_plddt_from_pdb(pdb2)

                        alignment  = align_sequences(seq1, seq2)
                        props      = compare_properties(
                            info or {}, info2 or {},
                            seq1, seq2, plddt1, plddt2
                        )
                        similarity = classify_similarity(alignment["identity"])

                    st.session_state.comparison_result = {
                        "seq1":      seq1,
                        "seq2":      seq2,
                        "plddt1":    plddt1,
                        "plddt2":    plddt2,
                        "alignment": alignment,
                        "props":     props,
                        "similarity": similarity,
                    }

        # ── Show comparison results ───────────────────────────
        if st.session_state.comparison_result:
            cr   = st.session_state.comparison_result
            aln  = cr["alignment"]
            props = cr["props"]
            sim  = cr["similarity"]
            info2 = st.session_state.protein_info2

            name1 = info["name"][:25] if info else st.session_state.loaded_id
            name2 = info2["name"][:25] if info2 else st.session_state.loaded_id2

            # ── Similarity banner ─────────────────────────────
            st.markdown(
                f"<div style='background:{sim['color']}22;"
                f"border:2px solid {sim['color']};"
                f"border-radius:12px; padding:20px; margin:16px 0; text-align:center'>"
                f"<h2 style='margin:0; color:{sim['color']}'>"
                f"{sim['emoji']} {aln['identity']}% Sequence Identity</h2>"
                f"<p style='color:{sim['color']}; margin:4px 0'>"
                f"<b>{sim['level']}</b></p>"
                f"<p style='color:#aaa; margin:8px 0 0 0'>{sim['description']}</p>"
                f"</div>",
                unsafe_allow_html=True
            )

            # ── Side-by-side 3D viewers ───────────────────────
            st.markdown("#### 🧬 3D Structures Side by Side")
            v1, v2 = st.columns(2)

            with v1:
                st.markdown(
                    f"<div style='text-align:center; color:#00d4ff; "
                    f"font-weight:bold; margin-bottom:8px'>A: {name1}</div>",
                    unsafe_allow_html=True
                )
                render_protein(pdb_data, style="cartoon", color="confidence")

            with v2:
                st.markdown(
                    f"<div style='text-align:center; color:#ff6b6b; "
                    f"font-weight:bold; margin-bottom:8px'>B: {name2}</div>",
                    unsafe_allow_html=True
                )
                render_protein(
                    st.session_state.pdb_data2,
                    style="cartoon",
                    color="confidence"
                )

            # ── Property comparison table ─────────────────────
            st.markdown("---")
            st.markdown("#### 📊 Property Comparison")

            p1, p2, p3 = st.columns(3)

            with p1:
                st.markdown(f"**🔵 {name1}**")
                st.metric("Length",       f"{props['length1']} aa")
                st.metric("Avg pLDDT",    f"{props['avg_plddt1']}")
                st.metric("High conf %",  f"{props['high_conf_pct1']}%")
                st.metric("Net charge",   f"{props['charge1']:+d}")
                st.metric("Avg hydrophobicity",
                          f"{props['hydrophobicity1']}")

            with p2:
                st.markdown("**⚖️ Difference**")
                len_diff = props["length1"] - props["length2"]
                plddt_diff = props["avg_plddt1"] - props["avg_plddt2"]
                charge_diff = props["charge1"] - props["charge2"]
                hydro_diff = props["hydrophobicity1"] - props["hydrophobicity2"]

                def diff_color(val):
                    return "#2ecc71" if abs(val) < 5 else "#e67e22" if abs(val) < 20 else "#e74c3c"

                st.markdown("<br>", unsafe_allow_html=True)
                for label, val, unit in [
                    ("Length",  len_diff,    " aa"),
                    ("pLDDT",   plddt_diff,  ""),
                    ("Conf %",  props["high_conf_pct1"] - props["high_conf_pct2"], "%"),
                    ("Charge",  charge_diff, ""),
                    ("Hydrophobicity", round(hydro_diff, 2), ""),
                ]:
                    color = "#2ecc71" if abs(val) < 2 else "#e67e22"
                    st.markdown(
                        f"<div style='text-align:center; color:{color}; "
                        f"font-size:0.9em; padding:6px 0'>"
                        f"{val:+.1f}{unit}</div>",
                        unsafe_allow_html=True
                    )

            with p3:
                st.markdown(f"**🔴 {name2}**")
                st.metric("Length",       f"{props['length2']} aa")
                st.metric("Avg pLDDT",    f"{props['avg_plddt2']}")
                st.metric("High conf %",  f"{props['high_conf_pct2']}%")
                st.metric("Net charge",   f"{props['charge2']:+d}")
                st.metric("Avg hydrophobicity",
                          f"{props['hydrophobicity2']}")

            # ── pLDDT overlay chart ───────────────────────────
            st.markdown("---")
            st.markdown("#### 📈 Confidence Score Comparison (pLDDT)")
            st.caption("Compare structural confidence across both proteins")

            fig_plddt = go.Figure()
            fig_plddt.add_trace(go.Scatter(
                y=cr["plddt1"],
                x=list(range(1, len(cr["plddt1"]) + 1)),
                name=f"A: {name1}",
                line=dict(color="#00d4ff", width=1.5),
                fill="tozeroy",
                fillcolor="rgba(0, 212, 255, 0.1)",
                hovertemplate="Residue %{x}<br>pLDDT: %{y:.1f}<extra>A</extra>"
            ))
            fig_plddt.add_trace(go.Scatter(
                y=cr["plddt2"],
                x=list(range(1, len(cr["plddt2"]) + 1)),
                name=f"B: {name2}",
                line=dict(color="#ff6b6b", width=1.5),
                fill="tozeroy",
                fillcolor="rgba(255, 107, 107, 0.1)",
                hovertemplate="Residue %{x}<br>pLDDT: %{y:.1f}<extra>B</extra>"
            ))
            fig_plddt.add_hline(
                y=70, line_dash="dot",
                line_color="#aaa",
                annotation_text="Confident threshold (70)"
            )
            fig_plddt.update_layout(
                height=320,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis_title="Residue Position",
                yaxis_title="pLDDT Score",
                yaxis_range=[0, 105],
                font=dict(color="white"),
                legend=dict(
                    bgcolor="rgba(0,0,0,0.5)",
                    bordercolor="#333"
                )
            )
            st.plotly_chart(fig_plddt, use_container_width=True)

            # ── Amino acid group distribution ─────────────────
            st.markdown("---")
            st.markdown("#### 🧪 Amino Acid Composition")

            groups  = ["nonpolar", "polar", "positive", "negative"]
            colors  = ["#3498db", "#2ecc71", "#e74c3c", "#e67e22"]
            vals1   = [props["groups1"].get(g, 0) for g in groups]
            vals2   = [props["groups2"].get(g, 0) for g in groups]
            labels  = ["Nonpolar", "Polar", "Positive", "Negative"]

            fig_comp = go.Figure()
            fig_comp.add_trace(go.Bar(
                name=f"A: {name1}",
                x=labels, y=vals1,
                marker_color="#00d4ff",
                text=[f"{v:.1f}%" for v in vals1],
                textposition="outside",
            ))
            fig_comp.add_trace(go.Bar(
                name=f"B: {name2}",
                x=labels, y=vals2,
                marker_color="#ff6b6b",
                text=[f"{v:.1f}%" for v in vals2],
                textposition="outside",
            ))
            fig_comp.update_layout(
                barmode="group",
                height=350,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                xaxis_title="Amino Acid Group",
                yaxis_title="Percentage (%)",
                font=dict(color="white"),
                legend=dict(bgcolor="rgba(0,0,0,0.5)"),
            )
            st.plotly_chart(fig_comp, use_container_width=True)

            # ── Sequence differences ──────────────────────────
            st.markdown("---")
            st.markdown("#### 🔎 Sequence Differences")
            st.markdown(
                f"Found **{aln['diff_count']}** differences out of "
                f"**{aln['length']}** aligned positions "
                f"({100 - aln['identity']:.1f}% divergence)"
            )

            if aln["differences"]:
                import pandas as pd
                diff_rows = []
                for d in aln["differences"][:50]:
                    aa1_props = AA_PROPERTIES.get(d["aa1"], {})
                    aa2_props = AA_PROPERTIES.get(d["aa2"], {})
                    charge_change = aa2_props.get("charge", 0) - aa1_props.get("charge", 0)
                    diff_rows.append({
                        "Pos A": d["pos1"],
                        "Pos B": d["pos2"],
                        f"AA in {name1[:10]}": d["aa1"],
                        f"AA in {name2[:10]}": d["aa2"],
                        "Charge Δ": f"{charge_change:+d}",
                        "Size Δ": aa2_props.get("size", 0) - aa1_props.get("size", 0),
                        "Group A": aa1_props.get("group", "?"),
                        "Group B": aa2_props.get("group", "?"),
                    })

                df_diff = pd.DataFrame(diff_rows)
                st.dataframe(df_diff, use_container_width=True, hide_index=True)

                if aln["diff_count"] > 50:
                    st.caption(f"Showing first 50 of {aln['diff_count']} differences")

                # Export differences
                csv_diff = df_diff.to_csv(index=False)
                st.download_button(
                    "📥 Download Differences CSV",
                    data=csv_diff,
                    file_name=f"comparison_{st.session_state.loaded_id}_vs_{st.session_state.loaded_id2}.csv",
                    mime="text/csv",
                )
            else:
                st.success("✅ Sequences are identical!")

            # ── Protein info comparison ───────────────────────
            st.markdown("---")
            st.markdown("#### 📋 Full Protein Info")
            info_c1, info_c2 = st.columns(2)

            with info_c1:
                st.markdown(f"**🔵 Protein A — {name1}**")
                if info:
                    st.markdown(f"**Organism:** {info.get('organism','?')}")
                    st.markdown(f"**Gene:** {info.get('gene_symbol','?')}")
                    st.markdown(f"**Length:** {info.get('length','?')} aa")
                    st.info(info.get("function","")[:300])
                st.markdown(
                    f"[UniProt ↗](https://www.uniprot.org/uniprot/{st.session_state.loaded_id})"
                )

            with info_c2:
                st.markdown(f"**🔴 Protein B — {name2}**")
                if info2:
                    st.markdown(f"**Organism:** {info2.get('organism','?')}")
                    st.markdown(f"**Gene:** {info2.get('gene_symbol','?')}")
                    st.markdown(f"**Length:** {info2.get('length','?')} aa")
                    st.info(info2.get("function","")[:300])
                st.markdown(
                    f"[UniProt ↗](https://www.uniprot.org/uniprot/{st.session_state.loaded_id2})"
                )

        else:
            st.markdown("""
            <div style='text-align:center; padding:40px; color:#aaa;
            background:#1a1a2e; border-radius:8px'>
                <h3>👆 Select Protein B above to compare</h3>
                <p>Try comparing <b>Alpha Hemoglobin (P69905)</b>
                vs <b>Beta Hemoglobin (P68871)</b> —
                two very similar proteins that differ at just a few positions</p>
            </div>
            """, unsafe_allow_html=True)
    # ════════════════════════════════════════════════════════════
    # TAB 5 — How to Use
    # ════════════════════════════════════════════════════════════
    with tab6:
        st.subheader("ℹ️ How to Use This App")

        st.markdown("""
        ## 🚀 Getting Started

        ### Step 1 — Load a Protein
        - Enter a **UniProt ID** in the sidebar (e.g. `P69905` for Hemoglobin)
        - Or click one of the quick-start buttons on the home page
        - The app fetches the 3D structure from **DeepMind's AlphaFold database**

        ### Step 2 — Explore the Structure (Structure tab)
        - Rotate and zoom the 3D viewer by clicking and dragging
        - Change the **color scheme** to see confidence scores (pLDDT)
        - Blue regions = high confidence | Red regions = flexible/uncertain
        - The chart below shows confidence per amino acid position

        ### Step 3 — Simulate a Mutation (Mutation Lab tab)
        - Enter a **residue position** (e.g. 6)
        - Select a **new amino acid** (e.g. Valine)
        - Click **Analyze Mutation**
        - The app calculates **ΔΔG** — how much the mutation destabilizes the protein
        - The mutated position glows red in the 3D viewer

        ### Step 4 — Predict Disease (Disease Prediction tab)
        - After analyzing a mutation, go to Disease Prediction
        - Click **Predict Disease**
        - The app runs 5 analyses simultaneously:
            1. 🤖 **XGBoost ML model** — trained on 376K real clinical variants
            2. 📋 **ClinVar** — NIH database of known disease-causing mutations
            3. 🧬 **SIFT** — evolutionary conservation score
            4. 🔬 **PolyPhen-2** — structural disruption score
            5. 🤖 **AlphaMissense** — DeepMind's AI pathogenicity model
        - A **Groq LLaMA AI** generates a plain-English explanation

        ---

        ## 🧬 Try These Famous Mutations

        | Protein | UniProt | Position | Change | Disease |
        |---------|---------|----------|--------|---------|
        | Hemoglobin | P69905 | 6 | → V (Valine) | Sickle Cell Anemia |
        | BRCA1 | P38398 | 61 | → G (Glycine) | Breast Cancer |
        | TP53 | P04637 | 175 | → H (Histidine) | Li-Fraumeni Syndrome |
        | KRAS | P01116 | 12 | → D (Aspartate) | Pancreatic Cancer |
        | LRRK2 | Q5S007 | 2019 | → S (Serine) | Parkinson's Disease |

        ---

        ## 🔬 Understanding the Scores

        | Score | Range | Meaning |
        |-------|-------|---------|
        | pLDDT | 0–100 | AlphaFold confidence in structure |
        | ΔΔG | kcal/mol | Stability change (positive = destabilizing) |
        | SIFT | 0–1 | 0 = deleterious, 1 = tolerated |
        | PolyPhen | 0–1 | 1 = probably damaging |
        | AlphaMissense | 0–1 | 1 = likely pathogenic |

        ---

        ## 📚 Data Sources
        - **AlphaFold** — DeepMind's protein structure database (200M structures)
        - **ClinVar** — NIH database of clinically significant variants
        - **Ensembl VEP** — EMBL-EBI variant effect predictor
        - **AlphaMissense** — DeepMind (2023) pathogenicity predictions
        - **DynaMut2** — University of Queensland stability predictor
        - **Groq LLaMA 3** — AI explanation generation
        """)