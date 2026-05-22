import py3Dmol
import streamlit.components.v1 as components

def render_protein(pdb_data: str, style: str = "cartoon", color: str = "confidence"):
    """Render an interactive 3D protein viewer using raw HTML."""

    view = py3Dmol.view(width=700, height=500)
    view.addModel(pdb_data, "pdb")

    if style == "cartoon":
        if color == "confidence":
            view.setStyle({"cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "roygb",
                    "min": 50,
                    "max": 100
                }
            }})
        elif color == "rainbow":
            view.setStyle({"cartoon": {"color": "spectrum"}})
        elif color == "secondary":
            view.setStyle({"cartoon": {"colorscheme": "ssPyMol"}})

    elif style == "stick":
        view.setStyle({"stick": {}})

    elif style == "sphere":
        view.setStyle({"sphere": {"colorscheme": "Jmol"}})

    elif style == "surface":
        view.setStyle({"cartoon": {"color": "white"}})
        view.addSurface(py3Dmol.VDW, {"opacity": 0.7, "colorscheme": "whiteCarbon"})

    view.setBackgroundColor("#0f0f0f")
    view.zoomTo()
    view.spin(True)

    # ← This is the key change: render as raw HTML instead of stmol
    html = view._make_html()
    components.html(html, height=520, scrolling=False)

def render_protein_with_mutation(pdb_data: str, mutation_position: int):
    """
    Render the protein in grey, with the mutation site
    highlighted as a glowing red sphere.
    """

    view = py3Dmol.view(width=700, height=500)
    view.addModel(pdb_data, "pdb")

    # Whole protein in grey cartoon
    view.setStyle({"cartoon": {"color": "#aaaaaa"}})

    # Mutation site — red sphere, much bigger than normal
    view.addStyle(
        {"resi": mutation_position},
        {"sphere": {"color": "#e74c3c", "radius": 1.2}}
    )

    # Also highlight the surrounding region (±3 residues) in orange
    for nearby in range(mutation_position - 3, mutation_position + 4):
        if nearby != mutation_position:
            view.addStyle(
                {"resi": nearby},
                {"cartoon": {"color": "#e67e22"}}
            )

    view.setBackgroundColor("#0f0f0f")
    view.zoomTo({"resi": mutation_position})   # zoom into mutation site
    view.spin(False)                            # stop spinning so user can inspect

    html = view._make_html()
    components.html(html, height=520, scrolling=False)