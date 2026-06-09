"""
protein_comparison.py
Handles all comparison logic between two proteins.
Sequence alignment, structural similarity, property differences.
"""

import requests
import urllib3
urllib3.disable_warnings()


# ── Amino acid properties for comparison ─────────────────────────
AA_PROPERTIES = {
    "A": {"hydrophobicity": 1.8,  "charge": 0,  "size": 1,
          "group": "nonpolar",   "name": "Alanine"},
    "C": {"hydrophobicity": 2.5,  "charge": 0,  "size": 2,
          "group": "polar",      "name": "Cysteine"},
    "D": {"hydrophobicity": -3.5, "charge": -1, "size": 2,
          "group": "negative",   "name": "Aspartate"},
    "E": {"hydrophobicity": -3.5, "charge": -1, "size": 3,
          "group": "negative",   "name": "Glutamate"},
    "F": {"hydrophobicity": 2.8,  "charge": 0,  "size": 4,
          "group": "nonpolar",   "name": "Phenylalanine"},
    "G": {"hydrophobicity": -0.4, "charge": 0,  "size": 1,
          "group": "nonpolar",   "name": "Glycine"},
    "H": {"hydrophobicity": -3.2, "charge": 1,  "size": 3,
          "group": "positive",   "name": "Histidine"},
    "I": {"hydrophobicity": 4.5,  "charge": 0,  "size": 3,
          "group": "nonpolar",   "name": "Isoleucine"},
    "K": {"hydrophobicity": -3.9, "charge": 1,  "size": 4,
          "group": "positive",   "name": "Lysine"},
    "L": {"hydrophobicity": 3.8,  "charge": 0,  "size": 3,
          "group": "nonpolar",   "name": "Leucine"},
    "M": {"hydrophobicity": 1.9,  "charge": 0,  "size": 3,
          "group": "nonpolar",   "name": "Methionine"},
    "N": {"hydrophobicity": -3.5, "charge": 0,  "size": 2,
          "group": "polar",      "name": "Asparagine"},
    "P": {"hydrophobicity": -1.6, "charge": 0,  "size": 2,
          "group": "nonpolar",   "name": "Proline"},
    "Q": {"hydrophobicity": -3.5, "charge": 0,  "size": 3,
          "group": "polar",      "name": "Glutamine"},
    "R": {"hydrophobicity": -4.5, "charge": 1,  "size": 5,
          "group": "positive",   "name": "Arginine"},
    "S": {"hydrophobicity": -0.8, "charge": 0,  "size": 2,
          "group": "polar",      "name": "Serine"},
    "T": {"hydrophobicity": -0.7, "charge": 0,  "size": 2,
          "group": "polar",      "name": "Threonine"},
    "V": {"hydrophobicity": 4.2,  "charge": 0,  "size": 2,
          "group": "nonpolar",   "name": "Valine"},
    "W": {"hydrophobicity": -0.9, "charge": 0,  "size": 5,
          "group": "nonpolar",   "name": "Tryptophan"},
    "Y": {"hydrophobicity": -1.3, "charge": 0,  "size": 4,
          "group": "polar",      "name": "Tyrosine"},
}


# ── Extract sequence from PDB ─────────────────────────────────────

def extract_sequence_from_pdb(pdb_data: str) -> str:
    """Extract the amino acid sequence from a PDB file."""
    three_to_one = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }

    sequence = []
    seen_residues = set()

    for line in pdb_data.split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num  = int(line[22:26].strip())
            res_name = line[17:20].strip()
            if res_num not in seen_residues:
                aa = three_to_one.get(res_name, "X")
                sequence.append(aa)
                seen_residues.add(res_num)

    return "".join(sequence)


def extract_plddt_from_pdb(pdb_data: str) -> list:
    """Extract pLDDT scores from PDB B-factor column."""
    scores = []
    seen   = set()
    for line in pdb_data.split("\n"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num = int(line[22:26].strip())
            if res_num not in seen:
                try:
                    scores.append(float(line[60:66].strip()))
                    seen.add(res_num)
                except:
                    pass
    return scores


# ── Sequence alignment ────────────────────────────────────────────

def align_sequences(seq1: str, seq2: str) -> dict:
    """
    Simple global sequence alignment using Needleman-Wunsch algorithm.
    Returns alignment strings, identity %, and difference positions.
    """
    n, m = len(seq1), len(seq2)

    # Scoring
    MATCH    =  1
    MISMATCH = -1
    GAP      = -2

    # Fill DP matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = i * GAP
    for j in range(m + 1):
        dp[0][j] = j * GAP

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match  = dp[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
            delete = dp[i-1][j]   + GAP
            insert = dp[i][j-1]   + GAP
            dp[i][j] = max(match, delete, insert)

    # Traceback
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = MATCH if seq1[i-1] == seq2[j-1] else MISMATCH
            if dp[i][j] == dp[i-1][j-1] + score:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif dp[i][j] == dp[i-1][j] + GAP:
                aligned1.append(seq1[i-1])
                aligned2.append("-")
                i -= 1
            else:
                aligned1.append("-")
                aligned2.append(seq2[j-1])
                j -= 1
        elif i > 0:
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            j -= 1

    a1 = "".join(reversed(aligned1))
    a2 = "".join(reversed(aligned2))

    # Calculate identity
    matches     = sum(1 for x, y in zip(a1, a2) if x == y and x != "-")
    aligned_len = sum(1 for x, y in zip(a1, a2) if x != "-" or y != "-")
    identity    = (matches / aligned_len * 100) if aligned_len > 0 else 0

    # Find difference positions (first 200 to keep it fast)
    differences = []
    pos1, pos2 = 0, 0
    for x, y in zip(a1, a2):
        if x != "-":
            pos1 += 1
        if y != "-":
            pos2 += 1
        if x != y and x != "-" and y != "-":
            differences.append({
                "pos1": pos1,
                "pos2": pos2,
                "aa1":  x,
                "aa2":  y,
            })
        if len(differences) >= 200:
            break

    return {
        "aligned1":   a1,
        "aligned2":   a2,
        "identity":   round(identity, 1),
        "matches":    matches,
        "length":     aligned_len,
        "differences": differences,
        "diff_count": len(differences),
    }


# ── Protein property comparison ───────────────────────────────────

def compare_properties(info1: dict, info2: dict,
                        seq1: str, seq2: str,
                        plddt1: list, plddt2: list) -> dict:
    """
    Compare two proteins across multiple dimensions.
    Returns a structured comparison dict.
    """

    # Basic properties
    len1 = len(seq1)
    len2 = len(seq2)

    avg_plddt1 = round(sum(plddt1) / len(plddt1), 1) if plddt1 else 0
    avg_plddt2 = round(sum(plddt2) / len(plddt2), 1) if plddt2 else 0

    high_conf1 = sum(1 for s in plddt1 if s >= 70) / len(plddt1) * 100 if plddt1 else 0
    high_conf2 = sum(1 for s in plddt2 if s >= 70) / len(plddt2) * 100 if plddt2 else 0

    # Amino acid composition
    def aa_composition(seq):
        comp = {}
        for aa in seq:
            comp[aa] = comp.get(aa, 0) + 1
        return {k: round(v / len(seq) * 100, 1) for k, v in comp.items()}

    comp1 = aa_composition(seq1)
    comp2 = aa_composition(seq2)

    # Chemical properties
    def calc_charge(seq):
        pos = sum(1 for aa in seq if AA_PROPERTIES.get(aa, {}).get("charge", 0) > 0)
        neg = sum(1 for aa in seq if AA_PROPERTIES.get(aa, {}).get("charge", 0) < 0)
        return pos - neg

    def calc_hydrophobicity(seq):
        total = sum(AA_PROPERTIES.get(aa, {}).get("hydrophobicity", 0) for aa in seq)
        return round(total / len(seq), 2) if seq else 0

    def calc_group_distribution(seq):
        groups = {"nonpolar": 0, "polar": 0, "positive": 0, "negative": 0}
        for aa in seq:
            g = AA_PROPERTIES.get(aa, {}).get("group", "nonpolar")
            groups[g] = groups.get(g, 0) + 1
        return {k: round(v / len(seq) * 100, 1) for k, v in groups.items()}

    return {
        # Basic
        "length1":         len1,
        "length2":         len2,
        "length_diff":     abs(len1 - len2),
        "longer":          info1.get("name", "Protein 1") if len1 > len2
                           else info2.get("name", "Protein 2"),

        # Confidence
        "avg_plddt1":      avg_plddt1,
        "avg_plddt2":      avg_plddt2,
        "high_conf_pct1":  round(high_conf1, 1),
        "high_conf_pct2":  round(high_conf2, 1),

        # Chemical
        "charge1":         calc_charge(seq1),
        "charge2":         calc_charge(seq2),
        "hydrophobicity1": calc_hydrophobicity(seq1),
        "hydrophobicity2": calc_hydrophobicity(seq2),

        # Group distribution
        "groups1":         calc_group_distribution(seq1),
        "groups2":         calc_group_distribution(seq2),

        # Composition
        "composition1":    comp1,
        "composition2":    comp2,
    }


# ── Similarity classification ─────────────────────────────────────

def classify_similarity(identity: float) -> dict:
    """
    Classify sequence identity into biological categories.
    Based on the twilight zone of sequence alignment.
    """
    if identity >= 90:
        return {
            "level":       "Nearly Identical",
            "color":       "#00d4ff",
            "emoji":       "🔵",
            "description": "These proteins are essentially the same — likely "
                           "the same protein from different sources or very "
                           "close variants.",
        }
    elif identity >= 70:
        return {
            "level":       "Highly Similar",
            "color":       "#2ecc71",
            "emoji":       "🟢",
            "description": "These proteins share high sequence similarity. "
                           "They likely have the same function and very "
                           "similar 3D structures.",
        }
    elif identity >= 40:
        return {
            "level":       "Moderately Similar",
            "color":       "#f1c40f",
            "emoji":       "🟡",
            "description": "These proteins share moderate similarity. "
                           "They may have related functions or be part "
                           "of the same protein family.",
        }
    elif identity >= 20:
        return {
            "level":       "Distantly Related",
            "color":       "#e67e22",
            "emoji":       "🟠",
            "description": "Low sequence similarity. These may be very "
                           "distantly related or share only a small "
                           "functional domain.",
        }
    else:
        return {
            "level":       "Unrelated",
            "color":       "#e74c3c",
            "emoji":       "🔴",
            "description": "These proteins appear to be unrelated based "
                           "on sequence alone. They may have completely "
                           "different structures and functions.",
        }


if __name__ == "__main__":
    print("Testing protein comparison...")
    seq1 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLG"
    seq2 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDD"
    result = align_sequences(seq1, seq2)
    print(f"Identity: {result['identity']}%")
    print(f"Differences: {result['diff_count']}")
    sim = classify_similarity(result['identity'])
    print(f"Classification: {sim['emoji']} {sim['level']}")