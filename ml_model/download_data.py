import pandas as pd
import os
import re

# ── Amino acid properties ─────────────────────────────────────────────────────
AA_PROPERTIES = {
    "A": {"hydrophobicity": 1.8,  "charge": 0,  "size": 1},
    "C": {"hydrophobicity": 2.5,  "charge": 0,  "size": 2},
    "D": {"hydrophobicity": -3.5, "charge": -1, "size": 2},
    "E": {"hydrophobicity": -3.5, "charge": -1, "size": 3},
    "F": {"hydrophobicity": 2.8,  "charge": 0,  "size": 4},
    "G": {"hydrophobicity": -0.4, "charge": 0,  "size": 1},
    "H": {"hydrophobicity": -3.2, "charge": 1,  "size": 3},
    "I": {"hydrophobicity": 4.5,  "charge": 0,  "size": 3},
    "K": {"hydrophobicity": -3.9, "charge": 1,  "size": 4},
    "L": {"hydrophobicity": 3.8,  "charge": 0,  "size": 3},
    "M": {"hydrophobicity": 1.9,  "charge": 0,  "size": 3},
    "N": {"hydrophobicity": -3.5, "charge": 0,  "size": 2},
    "P": {"hydrophobicity": -1.6, "charge": 0,  "size": 2},
    "Q": {"hydrophobicity": -3.5, "charge": 0,  "size": 3},
    "R": {"hydrophobicity": -4.5, "charge": 1,  "size": 5},
    "S": {"hydrophobicity": -0.8, "charge": 0,  "size": 2},
    "T": {"hydrophobicity": -0.7, "charge": 0,  "size": 2},
    "V": {"hydrophobicity": 4.2,  "charge": 0,  "size": 2},
    "W": {"hydrophobicity": -0.9, "charge": 0,  "size": 5},
    "Y": {"hydrophobicity": -1.3, "charge": 0,  "size": 4},
}


def load_local_clinvar(filepath: str) -> pd.DataFrame:
    """Load the manually downloaded variant_summary.txt file."""

    print(f"📂 Loading local file: {filepath}")

    if not os.path.exists(filepath):
        print(f"❌ File not found at: {filepath}")
        print("   Make sure variant_summary.txt is inside the ml_model/ folder")
        return None

    size_mb = os.path.getsize(filepath) / (1024 * 1024)
    print(f"   File size: {size_mb:.0f} MB")
    print("   Reading... (may take 30–60 seconds for a large file)")

    df = pd.read_csv(
        filepath,
        sep="\t",
        low_memory=False,
        on_bad_lines="skip"   # skip any malformed rows
    )

    print(f"✅ Loaded {len(df):,} variants")
    print(f"   Columns: {list(df.columns[:8])}...")  # show first 8 columns
    return df


def prepare_training_data(df: pd.DataFrame) -> pd.DataFrame:
    """Filter down to clean, labeled missense variants."""

    print("\n🔧 Preparing training data...")

    # Keep only single nucleotide variants (missense mutations)
    df = df[df["Type"] == "single nucleotide variant"].copy()
    print(f"   After filtering SNVs: {len(df):,}")

    # Keep only clear pathogenic / benign labels
    label_map = {
        "Pathogenic":        1,
        "Likely pathogenic": 1,
        "Benign":            0,
        "Likely benign":     0,
    }
    df = df[df["ClinicalSignificance"].isin(label_map.keys())].copy()
    df["label"] = df["ClinicalSignificance"].map(label_map)
    print(f"   After filtering clear labels: {len(df):,}")

    # Keep useful columns only
    keep_cols = [
        "GeneSymbol", "Name", "ClinicalSignificance",
        "label", "ReviewStatus", "NumberSubmitters",
        "Chromosome", "PhenotypeList", "Origin",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df[keep_cols].copy()

    # Clean up
    df = df.dropna(subset=["label"])
    df["NumberSubmitters"] = pd.to_numeric(
        df.get("NumberSubmitters", 1), errors="coerce"
    ).fillna(1)

    print(f"   Final dataset: {len(df):,} variants")
    print(f"\n   Label distribution:")
    print(f"   🔴 Pathogenic: {(df['label']==1).sum():,}")
    print(f"   🟢 Benign:     {(df['label']==0).sum():,}")

    return df


def extract_amino_acids_from_name(name: str):
    """
    Parse ClinVar mutation name to get original and new amino acid.
    Handles formats like:
      p.Glu7Val   (3-letter)
      p.E7V       (1-letter)
    """
    three_to_one = {
        "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F",
        "Gly": "G", "His": "H", "Ile": "I", "Lys": "K", "Leu": "L",
        "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q", "Arg": "R",
        "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
    }

    name = str(name)

    # 3-letter format: p.Glu7Val
    match = re.search(r'p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})', name)
    if match:
        orig = three_to_one.get(match.group(1))
        new  = three_to_one.get(match.group(2))
        if orig and new:
            return orig, new

    # 1-letter format: p.E7V
    match = re.search(r'p\.([A-Z])\d+([A-Z])', name)
    if match:
        return match.group(1), match.group(2)

    return None, None


def extract_features(df: pd.DataFrame) -> pd.DataFrame:
    """Convert each variant into numerical features for the ML model."""

    print("\n⚙️  Extracting features...")
    print("   (Processing rows — this takes a minute...)")

    features = []

    for i, (_, row) in enumerate(df.iterrows()):

        # Progress update every 50,000 rows
        if i % 50000 == 0 and i > 0:
            print(f"   Processed {i:,} rows...")

        orig_aa, new_aa = extract_amino_acids_from_name(row.get("Name", ""))

        if orig_aa in AA_PROPERTIES and new_aa in AA_PROPERTIES:
            orig = AA_PROPERTIES[orig_aa]
            new  = AA_PROPERTIES[new_aa]

            hydro_diff  = abs(orig["hydrophobicity"] - new["hydrophobicity"])
            charge_diff = abs(orig["charge"]         - new["charge"])
            size_diff   = abs(orig["size"]            - new["size"])
            charge_flip = 1 if orig["charge"] != new["charge"] else 0
            has_aa_info = 1

        else:
            # Couldn't parse amino acids — use zeros
            hydro_diff  = 0.0
            charge_diff = 0
            size_diff   = 0
            charge_flip = 0
            has_aa_info = 0

        features.append({
            "hydrophobicity_diff": hydro_diff,
            "charge_diff":         charge_diff,
            "size_diff":           size_diff,
            "charge_flip":         charge_flip,
            "has_aa_info":         has_aa_info,
            "num_submitters":      float(row.get("NumberSubmitters", 1)),
            "label":               int(row["label"]),
        })

    feature_df = pd.DataFrame(features)

    # Only keep rows where we successfully parsed amino acids
    parsed = feature_df[feature_df["has_aa_info"] == 1]
    print(f"\n   Total rows:          {len(feature_df):,}")
    print(f"   Successfully parsed: {len(parsed):,}")
    print(f"   Skipped (no AA):     {len(feature_df) - len(parsed):,}")

    return parsed


if __name__ == "__main__":

    # ── Point to your local file ──────────────────────────────────
    LOCAL_FILE = "variant_summary.txt"

    # Load it
    df_raw = load_local_clinvar(LOCAL_FILE)
    if df_raw is None:
        exit()

    # Clean it
    df_clean = prepare_training_data(df_raw)

    # Extract features
    df_features = extract_features(df_clean)

    # Save both outputs
    df_clean.to_csv("clinvar_clean.csv", index=False)
    df_features.to_csv("training_features.csv", index=False)

    print("\n✅ Done! Saved:")
    print("   clinvar_clean.csv      ← full cleaned dataset")
    print("   training_features.csv  ← features ready for training")
    print("\n▶️  Next: python ml_model/train_model.py")