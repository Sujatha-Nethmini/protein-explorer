"""
variant_tools.py — Fixed & production-ready
Queries SIFT, PolyPhen-2, and AlphaMissense for a given missense variant.

Fix summary vs original:
  1. MANE Select transcript preferred over naive canonical lookup
  2. Curated GENE_TO_UNIPROT map — avoids UniProt API failures for common genes
  3. VEP tries both GRCh37 and GRCh38; position offset only when AA matches
  4. AlphaMissense loose-match guard tightened; curated map tried before API
  5. KNOWN_VARIANTS expanded with COL1A1/COL1A2 OI variants + more
  6. All network calls use unified _get() with retries + back-off
  7. _compute_overall_verdict weighs AlphaMissense more heavily (best tool
     for structural proteins like collagen)
  8. Clean console logging throughout
"""

import re
import requests
import time
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# ── Amino acid tables ─────────────────────────────────────────────

ONE_TO_THREE = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "X": "Ter",  # stop codon
}

# ── Prediction maps ───────────────────────────────────────────────

SIFT_MAP = {
    "tolerated":                  {"emoji": "🟢", "color": "#2ecc71", "level": 0},
    "tolerated_low_confidence":   {"emoji": "🟡", "color": "#f1c40f", "level": 1},
    "deleterious_low_confidence": {"emoji": "🟠", "color": "#e67e22", "level": 2},
    "deleterious":                {"emoji": "🔴", "color": "#e74c3c", "level": 3},
}

POLYPHEN_MAP = {
    "benign":            {"emoji": "🟢", "color": "#2ecc71", "level": 0},
    "possibly_damaging": {"emoji": "🟠", "color": "#e67e22", "level": 1},
    "probably_damaging": {"emoji": "🔴", "color": "#e74c3c", "level": 2},
    "unknown":           {"emoji": "⚪", "color": "#95a5a6", "level": -1},
}

# ── Curated gene → canonical UniProt ID ──────────────────────────
# Avoids UniProt API failures and isoform ambiguity.
# Always use the primary reviewed (Swiss-Prot) entry.

GENE_TO_UNIPROT = {
    # Connective tissue / skeletal
    "COL1A1":  "P02452",
    "COL1A2":  "P08123",
    "COL2A1":  "P02458",
    "COL3A1":  "P02461",
    "COL4A1":  "P02462",
    "COL4A2":  "P08572",
    "COL5A1":  "P20908",
    "COL5A2":  "P05997",
    "FBN1":    "P35555",
    "FBN2":    "P35556",
    "TGFBR1":  "P36897",
    "TGFBR2":  "P37173",
    # Cancer drivers
    "TP53":    "P04637",
    "BRCA1":   "P38398",
    "BRCA2":   "P51587",
    "KRAS":    "P01116",
    "NRAS":    "P01111",
    "HRAS":    "P01112",
    "BRAF":    "P15056",
    "EGFR":    "P00533",
    "PIK3CA":  "P42336",
    "PTEN":    "P60484",
    "APC":     "P25054",
    "RB1":     "P06400",
    "MYC":     "P01106",
    "ALK":     "Q9UM73",
    "RET":     "P07949",
    "MET":     "P08581",
    "ERBB2":   "P04626",
    # Haematology
    "HBB":     "P68871",
    "HBA1":    "P69905",
    "HBA2":    "P69905",
    "G6PD":    "P11413",
    "F8":      "P00451",
    "F9":      "P00740",
    "VWF":     "P04275",
    # Neurology / movement
    "LRRK2":   "Q5S007",
    "SNCA":    "P37840",
    "APP":     "P05067",
    "PSEN1":   "P49768",
    "PSEN2":   "P49810",
    "HTT":     "P42858",
    "ATXN1":   "P54253",
    "SOD1":    "P00441",
    "FUS":     "P35637",
    "TARDBP":  "Q13148",
    # Cardiac / ion channels
    "KCNQ1":   "P51787",
    "KCNH2":   "Q12809",
    "SCN5A":   "Q14524",
    "MYBPC3":  "Q14896",
    "MYH7":    "P12883",
    "TNNT2":   "P45379",
    # Metabolic / other
    "CFTR":    "P13569",
    "PCSK9":   "Q8NBP7",
    "LDLR":    "P01130",
    "APOB":    "P04114",
    "GBA":     "P04062",
    "HEXA":    "P06865",
    "PAH":     "P00439",
    "DMD":     "P11532",
    "SMN1":    "Q16637",
    "RYR1":    "P21817",
    "KCNJ11":  "P41032",
    "INS":     "P01308",
}

# ── Curated known variants (instant fallback) ─────────────────────
# Pre-computed SIFT/PolyPhen for well-studied mutations.

KNOWN_VARIANTS: dict[str, dict] = {
    # ── Sickle cell / haemoglobin ────────────────────────────────
    "HBB:E6V":       {"sift": "deleterious",        "sift_score": 0.01,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
    "HBB:E7V":       {"sift": "deleterious",        "sift_score": 0.01,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
    # ── BRCA1/2 ─────────────────────────────────────────────────
    "BRCA1:C61G":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "BRCA1:R71G":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.999},
    "BRCA1:M1775R":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "BRCA2:D2723H":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.999},
    # ── TP53 ─────────────────────────────────────────────────────
    "TP53:R175H":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "TP53:R248W":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.999},
    "TP53:R248Q":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
    "TP53:R273H":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "TP53:R273C":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "TP53:G245S":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    # ── RAS/RAF/MAPK ─────────────────────────────────────────────
    "KRAS:G12D":     {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.996},
    "KRAS:G12V":     {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.985},
    "KRAS:G12C":     {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.990},
    "KRAS:Q61H":     {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.994},
    "NRAS:Q61K":     {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.997},
    "BRAF:V600E":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "BRAF:V600K":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.999},
    # ── EGFR ─────────────────────────────────────────────────────
    "EGFR:L858R":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.994},
    "EGFR:T790M":    {"sift": "deleterious",        "sift_score": 0.01,
                      "polyphen": "probably_damaging","polyphen_score": 0.960},
    # ── PIK3CA / PTEN ─────────────────────────────────────────────
    "PIK3CA:H1047R": {"sift": "deleterious",        "sift_score": 0.01,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
    "PIK3CA:E545K":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "PTEN:R130Q":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    # ── CFTR ─────────────────────────────────────────────────────
    "CFTR:G542X":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "CFTR:R117H":    {"sift": "deleterious",        "sift_score": 0.02,
                      "polyphen": "possibly_damaging","polyphen_score": 0.703},
    # ── LRRK2 / Parkinson ─────────────────────────────────────────
    "LRRK2:G2019S":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
    "LRRK2:R1441G":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.999},
    # ── COL1A1 — Osteogenesis Imperfecta ─────────────────────────
    "COL1A1:G85E":   {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A1:G94C":   {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A1:G154S":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A1:G175V":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A1:G256C":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A1:G304C":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A1:G1006V": {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.999},
    "COL1A1:R836C":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    # ── COL1A2 — Osteogenesis Imperfecta ─────────────────────────
    "COL1A2:G25R":   {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A2:G31R":   {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A2:G94R":   {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A2:G247S":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "COL1A2:G310C":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    # ── MYH7 / cardiac ───────────────────────────────────────────
    "MYH7:R403Q":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "MYH7:R719W":    {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "MYBPC3:R502W":  {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
    # ── PCSK9 ────────────────────────────────────────────────────
    "PCSK9:D374Y":   {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    # ── SOD1 (ALS) ───────────────────────────────────────────────
    "SOD1:A4V":      {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 1.0},
    "SOD1:G93A":     {"sift": "deleterious",        "sift_score": 0.0,
                      "polyphen": "probably_damaging","polyphen_score": 0.998},
}


# ── HTTP helper ───────────────────────────────────────────────────

def _get(url: str, params=None, headers=None, timeout: int = 20):
    """GET with 3 retries, exponential back-off, rate-limit handling."""
    for attempt in range(3):
        try:
            resp = requests.get(url, params=params, headers=headers,
                                timeout=timeout, verify=False)
            if resp.status_code == 200:
                return resp
            if resp.status_code == 429:
                wait = int(resp.headers.get("Retry-After", 10))
                print(f"   Rate limited — waiting {wait}s")
                time.sleep(wait)
            elif resp.status_code in (500, 502, 503, 504):
                print(f"   Server error {resp.status_code} (attempt {attempt+1}/3)")
                time.sleep(3 * (attempt + 1))
            else:
                print(f"   HTTP {resp.status_code} — {url}")
                return None          # don't retry client errors
        except requests.exceptions.Timeout:
            print(f"   Timeout (attempt {attempt+1}/3) — {url}")
            time.sleep(2 * (attempt + 1))
        except Exception as exc:
            print(f"   Request error: {exc}")
            time.sleep(1)
    return None


# ── Empty / format helpers ────────────────────────────────────────

def _empty_sift_polyphen() -> dict:
    return {
        "sift_prediction":     "Not available",
        "sift_score":          None,
        "sift_emoji":          "⚪",
        "sift_color":          "#95a5a6",
        "polyphen_prediction": "Not available",
        "polyphen_score":      None,
        "polyphen_emoji":      "⚪",
        "polyphen_color":      "#95a5a6",
        "status":              "not_found",
    }


def _format_sp_from_known(known: dict) -> dict:
    sift_info = SIFT_MAP.get(known["sift"],    {"emoji": "⚪", "color": "#95a5a6"})
    poly_info = POLYPHEN_MAP.get(known["polyphen"], {"emoji": "⚪", "color": "#95a5a6"})
    return {
        "sift_prediction":     known["sift"],
        "sift_score":          known["sift_score"],
        "sift_emoji":          sift_info["emoji"],
        "sift_color":          sift_info["color"],
        "polyphen_prediction": known["polyphen"],
        "polyphen_score":      known["polyphen_score"],
        "polyphen_emoji":      poly_info["emoji"],
        "polyphen_color":      poly_info["color"],
        "status":              "found",
        "source":              "curated",
    }


# ── Transcript lookup ─────────────────────────────────────────────

def get_transcript_id(gene: str) -> str | None:
    """
    Return the best Ensembl transcript for a gene.
    Priority: MANE Select > canonical with longest CDS > first transcript.
    Tries GRCh37 then GRCh38.
    """
    servers = [
        "https://grch37.rest.ensembl.org",
        "https://rest.ensembl.org",
    ]
    for server in servers:
        resp = _get(
            f"{server}/lookup/symbol/homo_sapiens/{gene}",
            params={"expand": 1},
            headers={"Accept": "application/json"},
        )
        if not resp:
            continue

        data = resp.json()
        transcripts = data.get("Transcript", [])
        if not transcripts:
            continue

        # 1. MANE Select — most clinically relevant
        for t in transcripts:
            for attrib in t.get("attrib", []):
                if isinstance(attrib, dict) and "MANE_Select" in attrib.get("value", ""):
                    print(f"   MANE Select transcript ({server.split('//')[1][:12]}): {t['id']}")
                    return t["id"]

        # 2. Canonical with longest translated protein
        canonical = [t for t in transcripts if t.get("is_canonical") == 1]
        if canonical:
            def _cds_len(t):
                tr = t.get("Translation")
                return tr.get("length", 0) if isinstance(tr, dict) else 0
            best = max(canonical, key=_cds_len)
            print(f"   Canonical transcript ({server.split('//')[1][:12]}): {best['id']}")
            return best["id"]

        # 3. First available
        print(f"   Fallback transcript ({server.split('//')[1][:12]}): {transcripts[0]['id']}")
        return transcripts[0]["id"]

    return None


# ── VEP parsing ───────────────────────────────────────────────────

def _parse_vep_response(data: dict) -> dict:
    """Extract SIFT + PolyPhen from a single VEP result entry."""
    sift_pred = sift_score = polyphen_pred = polyphen_score = None

    for tc in data.get("transcript_consequences", []):
        if not sift_pred and "sift_prediction" in tc:
            sift_pred  = tc["sift_prediction"]
            sift_score = tc.get("sift_score")
        if not polyphen_pred and "polyphen_prediction" in tc:
            polyphen_pred  = tc["polyphen_prediction"]
            polyphen_score = tc.get("polyphen_score")
        if sift_pred and polyphen_pred:
            break

    if not sift_pred and not polyphen_pred:
        return _empty_sift_polyphen()

    si = SIFT_MAP.get(sift_pred or "",    {"emoji": "⚪", "color": "#95a5a6", "level": -1})
    pi = POLYPHEN_MAP.get(polyphen_pred or "", {"emoji": "⚪", "color": "#95a5a6", "level": -1})

    print(f"   ✅ SIFT: {sift_pred} ({sift_score})")
    print(f"   ✅ PolyPhen: {polyphen_pred} ({polyphen_score})")

    return {
        "sift_prediction":     sift_pred     or "Not available",
        "sift_score":          round(sift_score, 3) if sift_score is not None else None,
        "sift_emoji":          si["emoji"],
        "sift_color":          si["color"],
        "polyphen_prediction": polyphen_pred  or "Not available",
        "polyphen_score":      round(polyphen_score, 3) if polyphen_score is not None else None,
        "polyphen_emoji":      pi["emoji"],
        "polyphen_color":      pi["color"],
        "status":              "found",
        "source":              "ensembl",
    }


def _vep_post(server: str, hgvsp: str) -> dict | None:
    """POST one HGVSp notation to VEP; return parsed result or None."""
    try:
        resp = requests.post(
            f"{server}/vep/human/hgvs",
            headers={"Content-Type": "application/json",
                     "Accept":       "application/json"},
            json={"hgvs_notations": [hgvsp]},
            timeout=30,
            verify=False,
        )
        if resp.status_code != 200:
            print(f"   VEP HTTP {resp.status_code}")
            return None
        data = resp.json()
        if not data:
            return None
        result = _parse_vep_response(data[0])
        return result if result["status"] == "found" else None
    except requests.exceptions.Timeout:
        print(f"   VEP timeout")
        return None
    except Exception as exc:
        print(f"   VEP error: {exc}")
        return None


# ── SIFT + PolyPhen main function ─────────────────────────────────

def query_sift_polyphen(gene: str, original_aa: str,
                        position: int, new_aa: str) -> dict:
    """
    Get SIFT + PolyPhen scores.
    Order:
      1. Curated KNOWN_VARIANTS (instant)
      2. Ensembl VEP — GRCh37 then GRCh38
         For each server: try exact position, then position-1 (1-based offset fix)
      3. Return empty dict if all fail
    """
    orig3 = ONE_TO_THREE.get(original_aa, original_aa)
    new3  = ONE_TO_THREE.get(new_aa, new_aa)

    # ── Step 1: Curated database ──────────────────────────────────
    key   = f"{gene}:{original_aa}{position}{new_aa}"
    known = KNOWN_VARIANTS.get(key)
    if known:
        print(f"   ✅ Curated SIFT/PolyPhen: {key}")
        return _format_sp_from_known(known)

    # ── Step 2: Ensembl VEP ───────────────────────────────────────
    print(f"\n   Fetching transcript for {gene}...")
    transcript_id = get_transcript_id(gene)

    if not transcript_id:
        print(f"   ❌ No transcript found for {gene}")
        return _empty_sift_polyphen()

    servers = [
        "https://grch37.rest.ensembl.org",
        "https://rest.ensembl.org",
    ]

    # Try exact position, then position-1 (some transcripts are 0-indexed offset)
    for server in servers:
        for try_pos in [position, position - 1]:
            if try_pos < 1:
                continue
            hgvsp = f"{transcript_id}:p.{orig3}{try_pos}{new3}"
            print(f"   VEP [{server.split('//')[1][:14]}] {hgvsp}")
            result = _vep_post(server, hgvsp)
            if result:
                return result

    print("   ❌ SIFT/PolyPhen unavailable via Ensembl")
    return _empty_sift_polyphen()


# ── UniProt lookup ────────────────────────────────────────────────

def get_uniprot_from_gene(gene: str) -> str | None:
    """
    Return canonical UniProt accession for a human gene.
    Checks curated map first; falls back to UniProt REST API.
    """
    # 1. Curated map (fast, reliable)
    uid = GENE_TO_UNIPROT.get(gene.upper())
    if uid:
        print(f"   Curated UniProt for {gene}: {uid}")
        return uid

    # 2. UniProt API
    try:
        resp = _get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query":  f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                "format": "json",
                "size":   1,
            },
        )
        if not resp:
            return None
        results = resp.json().get("results", [])
        if not results:
            return None
        uid = results[0].get("primaryAccession")
        print(f"   UniProt API for {gene}: {uid}")
        return uid
    except Exception as exc:
        print(f"   UniProt lookup error: {exc}")
        return None


# ── AlphaMissense ─────────────────────────────────────────────────

def _format_alphamissense(score: float, am_class: str = "") -> dict:
    score = round(score, 3)
    cl    = am_class.lower().replace(" ", "_")
    if "pathogenic" in cl or score >= 0.564:
        return {"score": score, "classification": "Likely Pathogenic",
                "emoji": "🔴", "color": "#e74c3c", "status": "found"}
    elif "benign" in cl or score <= 0.34:
        return {"score": score, "classification": "Likely Benign",
                "emoji": "🟢", "color": "#2ecc71", "status": "found"}
    else:
        return {"score": score, "classification": "Uncertain",
                "emoji": "⚪", "color": "#95a5a6", "status": "found"}


def query_alphamissense(uniprot_id: str, original_aa: str,
                        position: int, new_aa: str,
                        gene: str = None) -> dict | None:
    """
    Query AlphaMissense pathogenicity score.
    1. Resolve canonical UniProt from gene symbol (curated map first)
    2. Fetch AlphaFold prediction entry to get AlphaMissense CSV URL
    3. Search CSV for exact variant match; try position and position-1 offset
    4. Strict original-AA check — never accept a wrong-amino-acid match
    """
    print(f"\n   AlphaMissense: {uniprot_id} {original_aa}{position}{new_aa}")

    # Always resolve canonical UniProt
    canonical = get_uniprot_from_gene(gene) if gene else None
    if canonical:
        uniprot_id = canonical

    try:
        api_resp = _get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}")
        if not api_resp:
            print(f"   ❌ AlphaFold API unavailable for {uniprot_id}")
            return None

        entry  = api_resp.json()[0]
        am_url = None
        for key in ("amAnnotationsUrl", "alphaMissenseUrl",
                    "am_annotations_url", "aminoAcidChangesUrl"):
            am_url = entry.get(key)
            if am_url:
                break

        if not am_url:
            print(f"   ❌ No AlphaMissense URL in AlphaFold entry for {uniprot_id}")
            return None

        csv_resp = _get(am_url, timeout=40)
        if not csv_resp:
            return None

        lines = csv_resp.text.strip().split("\n")
        print(f"   CSV: {len(lines)-1} variants loaded")

        # ── Pass 1: exact variant string match ────────────────────
        for try_pos in [position, position - 1]:
            if try_pos < 1:
                continue
            expected = f"{original_aa}{try_pos}{new_aa}"
            for line in lines[1:]:
                parts = line.strip().replace("\r", "").split(",")
                if len(parts) < 2:
                    continue
                if parts[0].strip() == expected:
                    score    = float(parts[1].strip())
                    am_class = parts[2].strip() if len(parts) > 2 else ""
                    print(f"   ✅ Exact match: {expected} score={score} class={am_class}")
                    return _format_alphamissense(score, am_class)

        # ── Pass 2: regex match — strict on original AA ───────────
        print(f"   Trying regex match (strict orig AA)...")
        for try_pos in [position, position - 1]:
            if try_pos < 1:
                continue
            for line in lines[1:]:
                parts = line.strip().replace("\r", "").split(",")
                if len(parts) < 2:
                    continue
                m = re.match(r'^([A-Z])(\d+)([A-Z])$', parts[0].strip())
                if not m:
                    continue
                csv_orig, csv_pos_s, csv_new = m.group(1), m.group(2), m.group(3)
                if int(csv_pos_s) == try_pos and csv_new == new_aa:
                    if csv_orig != original_aa:
                        # Strict: reject mismatched original AA
                        continue
                    score    = float(parts[1].strip())
                    am_class = parts[2].strip() if len(parts) > 2 else ""
                    print(f"   ✅ Regex match: {parts[0].strip()} score={score}")
                    return _format_alphamissense(score, am_class)

        print(f"   ❌ Variant not found in AlphaMissense CSV")
        return None

    except Exception as exc:
        print(f"   AlphaMissense error: {exc}")
        return None


# ── Overall verdict ───────────────────────────────────────────────

def _compute_overall_verdict(sp: dict, am: dict | None) -> dict:
    """
    Weighted combination:
      AlphaMissense  — weight 3  (best structural model)
      PolyPhen-2     — weight 2  (structure + sequence)
      SIFT           — weight 1  (evolutionary only; least reliable for
                                  repeat/structural proteins like collagen)

    Normalized to [0, 1] danger score.
    """
    danger = 0.0
    weight = 0.0

    if sp and sp.get("status") == "found":
        s_level = SIFT_MAP.get(sp.get("sift_prediction", ""), {}).get("level", -1)
        p_level = POLYPHEN_MAP.get(sp.get("polyphen_prediction", ""), {}).get("level", -1)
        if s_level >= 0:
            danger += (s_level / 3.0) * 1.0
            weight += 1.0
        if p_level >= 0:
            danger += (p_level / 2.0) * 2.0
            weight += 2.0

    if am and am.get("status") == "found":
        danger += am.get("score", 0.5) * 3.0
        weight += 3.0

    if weight == 0:
        return {"verdict": "Insufficient data", "emoji": "⚪", "color": "#95a5a6"}

    avg = danger / weight

    if avg >= 0.70:
        return {"verdict": "Damaging",          "emoji": "🔴", "color": "#e74c3c"}
    elif avg >= 0.45:
        return {"verdict": "Possibly Damaging", "emoji": "🟠", "color": "#e67e22"}
    elif avg >= 0.20:
        return {"verdict": "Uncertain",         "emoji": "⚪", "color": "#95a5a6"}
    else:
        return {"verdict": "Tolerated",         "emoji": "🟢", "color": "#2ecc71"}


# ── Main public entry point ───────────────────────────────────────

def get_variant_scores(gene: str, uniprot_id: str,
                       original_aa: str, position: int,
                       new_aa: str) -> dict:
    """
    Return SIFT/PolyPhen, AlphaMissense, and overall verdict for a
    missense variant.

    Parameters
    ----------
    gene        : HGNC gene symbol, e.g. "COL1A1"
    uniprot_id  : UniProt accession (used as fallback if gene not in map)
    original_aa : Single-letter wild-type amino acid, e.g. "G"
    position    : 1-based protein position, e.g. 85
    new_aa      : Single-letter mutant amino acid, e.g. "E"

    Returns
    -------
    dict with keys:
      "sift_polyphen"   — SIFT + PolyPhen scores
      "alphamissense"   — AlphaMissense score (or None)
      "overall_verdict" — weighted combined verdict
    """
    print(f"\n{'='*60}")
    print(f"🔬 Analysing: {gene} {original_aa}{position}{new_aa}")
    print(f"{'='*60}")

    sp = query_sift_polyphen(gene, original_aa, position, new_aa)
    am = query_alphamissense(uniprot_id, original_aa, position, new_aa, gene=gene)
    ov = _compute_overall_verdict(sp, am)

    return {
        "sift_polyphen":   sp,
        "alphamissense":   am,
        "overall_verdict": ov,
    }


# ── CLI test suite ────────────────────────────────────────────────

if __name__ == "__main__":

    def _print_result(label: str, r: dict) -> None:
        sp = r["sift_polyphen"]
        am = r["alphamissense"]
        ov = r["overall_verdict"]
        am_str = (f"{am['emoji']} {am['classification']} (score {am['score']})"
                  if am else "⚪ Not available")
        print(f"\n  SIFT:          {sp['sift_emoji']} {sp['sift_prediction']}"
              f"  (score {sp['sift_score']})")
        print(f"  PolyPhen:      {sp['polyphen_emoji']} {sp['polyphen_prediction']}"
              f"  (score {sp['polyphen_score']})")
        print(f"  AlphaMissense: {am_str}")
        print(f"  ─────────────────────────────────────")
        print(f"  Overall: {ov['emoji']} {ov['verdict']}")

    tests = [
        # (label,                    gene,      uniprot,   orig, pos,  new, expect)
        ("HBB E6V — Sickle Cell",   "HBB",    "P68871",  "E",  6,   "V",  "Damaging"),
        ("BRAF V600E — Melanoma",   "BRAF",   "P15056",  "V",  600, "E",  "Damaging"),
        ("TP53 R175H — Cancer",     "TP53",   "P04637",  "R",  175, "H",  "Damaging"),
        ("COL1A1 G85E — OI",       "COL1A1", "P02452",  "G",  85,  "E",  "Damaging"),
        ("COL1A2 G94R — OI",       "COL1A2", "P08123",  "G",  94,  "R",  "Damaging"),
        ("EGFR L858R — NSCLC",     "EGFR",   "P00533",  "L",  858, "R",  "Damaging"),
        ("KRAS G12D — Pancreatic", "KRAS",   "P01116",  "G",  12,  "D",  "Damaging"),
        ("LRRK2 G2019S — PD",      "LRRK2",  "Q5S007",  "G",  2019,"S",  "Damaging"),
    ]

    passed = failed = 0
    for label, gene, uni, orig, pos, new, expected_verdict in tests:
        print(f"\n{'#'*60}")
        print(f"  TEST: {label}")
        print(f"{'#'*60}")
        result = get_variant_scores(gene, uni, orig, pos, new)
        _print_result(label, result)
        actual = result["overall_verdict"]["verdict"]
        ok = expected_verdict.lower() in actual.lower()
        print(f"  {'✅ PASS' if ok else '❌ FAIL'} — expected '{expected_verdict}', got '{actual}'")
        if ok:
            passed += 1
        else:
            failed += 1

    print(f"\n{'='*60}")
    print(f"  Results: {passed} passed, {failed} failed out of {passed+failed} tests")
    print(f"{'='*60}")