import requests

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


def fetch_pdb(uniprot_id: str) -> str | None:
    """Get the correct PDB URL from AlphaFold's API, then download it."""

    # Step 1: Ask AlphaFold API for the correct download link
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    api_response = requests.get(api_url)

    if api_response.status_code != 200:
        return None

    data = api_response.json()

    if not data or len(data) == 0:
        return None

    # Step 2: Extract the actual PDB download URL from the response
    pdb_url = data[0].get("pdbUrl")

    if not pdb_url:
        return None

    # Step 3: Download the PDB file
    pdb_response = requests.get(pdb_url)

    if pdb_response.status_code == 200:
        return pdb_response.text

    return None


def fetch_protein_info(uniprot_id: str) -> dict | None:
    """Get protein name, organism, length, and function from UniProt."""
    url = f"{UNIPROT_BASE}/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code != 200:
        return None

    data = response.json()

    name = data.get("proteinDescription", {}) \
        .get("recommendedName", {}) \
        .get("fullName", {}) \
        .get("value", "Unknown")

    organism = data.get("organism", {}) \
        .get("scientificName", "Unknown")

    length = data.get("sequence", {}).get("length", 0)

    comments = data.get("comments", [])
    function = "No description available."
    for comment in comments:
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function = texts[0].get("value", function)
                break

    return {
        "name": name,
        "organism": organism,
        "length": length,
        "function": function,
        "id": uniprot_id
    }