import requests
import time

# ── Amino acid properties (kept as fallback) ─────────────────────────────────
AA_PROPERTIES = {
    "A": {"hydrophobicity": 1.8,  "charge": 0,  "size": 1, "name": "Alanine"},
    "C": {"hydrophobicity": 2.5,  "charge": 0,  "size": 2, "name": "Cysteine"},
    "D": {"hydrophobicity": -3.5, "charge": -1, "size": 2, "name": "Aspartate"},
    "E": {"hydrophobicity": -3.5, "charge": -1, "size": 3, "name": "Glutamate"},
    "F": {"hydrophobicity": 2.8,  "charge": 0,  "size": 4, "name": "Phenylalanine"},
    "G": {"hydrophobicity": -0.4, "charge": 0,  "size": 1, "name": "Glycine"},
    "H": {"hydrophobicity": -3.2, "charge": 1,  "size": 3, "name": "Histidine"},
    "I": {"hydrophobicity": 4.5,  "charge": 0,  "size": 3, "name": "Isoleucine"},
    "K": {"hydrophobicity": -3.9, "charge": 1,  "size": 4, "name": "Lysine"},
    "L": {"hydrophobicity": 3.8,  "charge": 0,  "size": 3, "name": "Leucine"},
    "M": {"hydrophobicity": 1.9,  "charge": 0,  "size": 3, "name": "Methionine"},
    "N": {"hydrophobicity": -3.5, "charge": 0,  "size": 2, "name": "Asparagine"},
    "P": {"hydrophobicity": -1.6, "charge": 0,  "size": 2, "name": "Proline"},
    "Q": {"hydrophobicity": -3.5, "charge": 0,  "size": 3, "name": "Glutamine"},
    "R": {"hydrophobicity": -4.5, "charge": 1,  "size": 5, "name": "Arginine"},
    "S": {"hydrophobicity": -0.8, "charge": 0,  "size": 2, "name": "Serine"},
    "T": {"hydrophobicity": -0.7, "charge": 0,  "size": 2, "name": "Threonine"},
    "V": {"hydrophobicity": 4.2,  "charge": 0,  "size": 2, "name": "Valine"},
    "W": {"hydrophobicity": -0.9, "charge": 0,  "size": 5, "name": "Tryptophan"},
    "Y": {"hydrophobicity": -1.3, "charge": 0,  "size": 4, "name": "Tyrosine"},
}

DYNAMUT2_URL = "http://biosig.lab.uq.edu.au/dynamut2/api/"


# ── Helpers ───────────────────────────────────────────────────────────────────

def get_residue_at_position(pdb_data: str, position: int) -> str | None:
    """Read the original amino acid letter at a given position."""
    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }
    for line in pdb_data.split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num = int(line[22:26].strip())
            if res_num == position:
                return three_to_one.get(line[17:20].strip())
    return None


def get_plddt_at_position(pdb_data: str, position: int) -> float | None:
    """Get the pLDDT confidence score at a given position."""
    for line in pdb_data.split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            if int(line[22:26].strip()) == position:
                try:
                    return float(line[60:66].strip())
                except:
                    return None
    return None


# ── DynaMut2 API (Approach A) ─────────────────────────────────────────────────

def query_dynamut2(pdb_data: str, chain: str, original_aa: str,
                   position: int, new_aa: str) -> dict | None:
    """
    Submit a mutation to the DynaMut2 server and poll for results.

    DynaMut2 is a proper molecular dynamics tool — it simulates
    actual atomic interactions to calculate ΔΔG, far more accurate
    than our local property-based estimate.

    Returns a dict with ddg and other scores, or None if it fails.
    """

    # Step 1 — Submit the job
    # Mutation format DynaMut2 expects: "CA6V"
    # = Chain + OriginalAA + Position + NewAA
    mutation_str = f"{chain}{original_aa}{position}{new_aa}"

    try:
        submit = requests.post(
            DYNAMUT2_URL,
            data={"mutation": mutation_str},
            files={"pdb_file": ("protein.pdb", pdb_data, "chemical/x-pdb")},
            timeout=30
        )

        if submit.status_code != 200:
            return None

        job_data = submit.json()

        # Step 2 — Poll for results (DynaMut2 is async — it processes in background)
        # It usually takes 10–30 seconds
        job_id = job_data.get("job_id") or job_data.get("id")

        if not job_id:
            # Some versions return results immediately
            return _parse_dynamut2_result(job_data)

        # Poll every 5 seconds, up to 60 seconds total
        for attempt in range(12):
            time.sleep(5)
            poll = requests.get(
                f"{DYNAMUT2_URL}results/{job_id}",
                timeout=15
            )
            if poll.status_code == 200:
                result = poll.json()
                if result.get("status") == "done" or "ddg" in result:
                    return _parse_dynamut2_result(result)

        return None  # timed out

    except Exception as e:
        print(f"DynaMut2 error: {e}")
        return None


def _parse_dynamut2_result(data: dict) -> dict | None:
    """Extract ΔΔG and related scores from DynaMut2 response."""
    try:
        # DynaMut2 returns multiple ΔΔG predictions — we average them
        ddg_values = []

        # Different keys depending on DynaMut2 version
        for key in ["ddg", "mcsm_ddg", "dynamut_ddg", "ddg_stability"]:
            val = data.get(key)
            if val is not None:
                try:
                    ddg_values.append(float(val))
                except:
                    pass

        if not ddg_values:
            return None

        avg_ddg = round(sum(ddg_values) / len(ddg_values), 2)

        return {
            "ddg":    avg_ddg,
            "source": "DynaMut2",
            "raw":    data,
            # DynaMut2 also gives us these bonus scores if available
            "solvation_ddg": data.get("solvation_ddg"),
            "vibrational_ddg": data.get("vibrational_ddg"),
        }

    except Exception as e:
        print(f"Parse error: {e}")
        return None


# ── Local fallback (Approach B) ───────────────────────────────────────────────

def calculate_ddg_local(original_aa: str, new_aa: str,
                        plddt: float = 80.0) -> dict:
    """
    Fallback: estimate ΔΔG from amino acid property differences.
    Used when DynaMut2 is unavailable.
    """
    orig = AA_PROPERTIES[original_aa]
    new  = AA_PROPERTIES[new_aa]

    hydrophobicity_diff = abs(orig["hydrophobicity"] - new["hydrophobicity"])
    charge_diff         = abs(orig["charge"]         - new["charge"])
    size_diff           = abs(orig["size"]            - new["size"])
    plddt_factor        = plddt / 100.0

    ddg = (
        (hydrophobicity_diff * 0.35) +
        (charge_diff         * 1.50) +
        (size_diff           * 0.40)
    ) * plddt_factor

    return {
        "ddg":                round(ddg, 2),
        "source":             "Local Estimate",
        "hydrophobicity_diff": round(hydrophobicity_diff, 2),
        "charge_diff":         charge_diff,
        "size_diff":           size_diff,
        "plddt_factor":        round(plddt_factor, 2),
        "original_name":       orig["name"],
        "new_name":            new["name"],
    }


# ── Main entry point ──────────────────────────────────────────────────────────

def calculate_ddg(pdb_data: str, original_aa: str, new_aa: str,
                  position: int, plddt: float = 80.0,
                  use_dynamut2: bool = True) -> dict:
    """
    Try DynaMut2 first. Fall back to local if it fails.
    Returns ddg_data dict with a 'source' key telling you which was used.
    """

    if use_dynamut2:
        result = query_dynamut2(pdb_data, "A", original_aa, position, new_aa)
        if result:
            # Add property info for the explanation layer
            orig = AA_PROPERTIES.get(original_aa, {})
            new  = AA_PROPERTIES.get(new_aa, {})
            result["hydrophobicity_diff"] = round(
                abs(orig.get("hydrophobicity", 0) - new.get("hydrophobicity", 0)), 2)
            result["charge_diff"]  = abs(orig.get("charge", 0) - new.get("charge", 0))
            result["size_diff"]    = abs(orig.get("size", 0)   - new.get("size", 0))
            result["original_name"] = orig.get("name", original_aa)
            result["new_name"]      = new.get("name", new_aa)
            return result

    # DynaMut2 failed or was skipped — use local
    return calculate_ddg_local(original_aa, new_aa, plddt)


# ── Interpretation (same for both approaches) ─────────────────────────────────

def interpret_ddg(ddg: float, original_aa: str, new_aa: str,
                  position: int) -> dict:
    """Turn ΔΔG into a human-readable result."""

    orig = AA_PROPERTIES[original_aa]
    new  = AA_PROPERTIES[new_aa]

    if ddg <= 0.3:
        level, color, emoji = "Neutral",                  "#2ecc71", "🟢"
    elif ddg <= 0.8:
        level, color, emoji = "Mildly Destabilizing",     "#f1c40f", "🟡"
    elif ddg <= 1.5:
        level, color, emoji = "Moderately Destabilizing", "#e67e22", "🟠"
    else:
        level, color, emoji = "Highly Destabilizing",     "#e74c3c", "🔴"

    summary = {
        "Neutral":                  "This mutation is likely harmless.",
        "Mildly Destabilizing":     "This mutation causes minor structural disruption.",
        "Moderately Destabilizing": "This mutation noticeably weakens the protein.",
        "Highly Destabilizing":     "This mutation severely disrupts protein stability.",
    }[level]

    parts = [
        f"At position {position}, the original amino acid is "
        f"**{orig['name']} ({original_aa})**."
    ]

    if orig["charge"] != new["charge"]:
        if orig["charge"] != 0 and new["charge"] == 0:
            parts.append(
                f"**{orig['name']}** carries an electrical charge that forms "
                f"important bonds. **{new['name']} ({new_aa})** is uncharged, "
                f"so those bonds break."
            )
        elif orig["charge"] == 0 and new["charge"] != 0:
            parts.append(
                f"Introducing a charged amino acid (**{new['name']}**) where "
                f"there was none can create electrostatic clashes."
            )
        else:
            parts.append(
                f"The charge flips, repelling or attracting nearby atoms "
                f"in unexpected ways."
            )

    if abs(orig["hydrophobicity"] - new["hydrophobicity"]) > 2.0:
        if orig["hydrophobicity"] < 0 and new["hydrophobicity"] > 0:
            parts.append(
                f"**{orig['name']}** is water-loving and sits on the surface. "
                f"**{new['name']}** is water-repelling and tries to hide inside "
                f"— creating a structural conflict."
            )
        elif orig["hydrophobicity"] > 0 and new["hydrophobicity"] < 0:
            parts.append(
                f"**{orig['name']}** is buried in the protein core. "
                f"Replacing it with water-loving **{new['name']}** "
                f"destabilizes that core."
            )

    if abs(orig["size"] - new["size"]) >= 2:
        if new["size"] > orig["size"]:
            parts.append(
                f"**{new['name']}** is significantly larger, causing steric "
                f"clashes — like forcing a big piece into a small puzzle slot."
            )
        else:
            parts.append(
                f"**{new['name']}** is much smaller, creating a void in the "
                f"protein core that destabilizes it."
            )

    parts.append(summary)

    return {
        "level":       level,
        "color":       color,
        "emoji":       emoji,
        "summary":     summary,
        "explanation": "\n\n".join(parts),
        "ddg":         ddg,
    }