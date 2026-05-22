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