import os
import ssl
import re
import time
import requests
import certifi
import urllib3
import urllib.request

# Fix macOS SSL
ssl._create_default_https_context = ssl.create_default_context
os.environ["SSL_CERT_FILE"] = certifi.where()
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

CLINVAR_SEARCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _get(url, params=None, timeout=30, retries=3):
    """Requests wrapper with retry logic for slow servers."""
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=timeout, verify=False)
            if resp.status_code == 200:
                return resp
            print(f"   HTTP {resp.status_code} (attempt {attempt+1}/{retries})")
        except requests.exceptions.Timeout:
            print(f"   Timeout (attempt {attempt+1}/{retries}) — retrying...")
            time.sleep(2 * (attempt + 1))  # wait longer each retry
        except Exception as e:
            print(f"   Request error: {e}")
            time.sleep(1)
    return None


def _fetch_variants_by_ids(id_list: list) -> dict:
    """Fetch ClinVar summaries — fetch in small batches to avoid timeouts."""
    time.sleep(0.5)

    all_results = {}

    # Fetch in batches of 3 instead of all at once
    for i in range(0, len(id_list), 3):
        batch = id_list[i:i+3]
        resp  = _get(CLINVAR_SUMMARY_URL, params={
            "db":      "clinvar",
            "id":      ",".join(batch),
            "retmode": "json",
        }, timeout=30)

        if resp and resp.status_code == 200:
            result = resp.json().get("result", {})
            all_results.update(result)

        time.sleep(0.5)  # be polite between batches

    return all_results


def _search_ids(term: str, retmax: int = 10) -> list:
    """Search ClinVar and return list of variant IDs."""
    resp = _get(CLINVAR_SEARCH_URL, params={
        "db":      "clinvar",
        "term":    term,
        "retmode": "json",
        "retmax":  retmax,
    })
    if not resp or resp.status_code != 200:
        return []
    data = resp.json()
    ids  = data.get("esearchresult", {}).get("idlist", [])
    cnt  = data.get("esearchresult", {}).get("count", 0)
    print(f"   '{term}' → {cnt} results")
    return ids


# ── Main ClinVar query ────────────────────────────────────────────────────────

ONE_TO_THREE = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
}

SIGNIFICANCE_MAP = {
    "Pathogenic":                    {"level": "pathogenic",        "emoji": "🔴", "color": "#e74c3c"},
    "Likely pathogenic":             {"level": "likely_pathogenic", "emoji": "🟠", "color": "#e67e22"},
    "Pathogenic/Likely pathogenic":  {"level": "likely_pathogenic", "emoji": "🟠", "color": "#e67e22"},
    "Benign":                        {"level": "benign",            "emoji": "🟢", "color": "#2ecc71"},
    "Likely benign":                 {"level": "likely_benign",     "emoji": "🟡", "color": "#f1c40f"},
    "Benign/Likely benign":          {"level": "likely_benign",     "emoji": "🟡", "color": "#f1c40f"},
    "Uncertain significance":        {"level": "uncertain",         "emoji": "⚪", "color": "#95a5a6"},
    "Conflicting interpretations":   {"level": "conflicting",       "emoji": "🔵", "color": "#3498db"},
}

PRIORITY = {
    "pathogenic": 0, "likely pathogenic": 1,
    "pathogenic/likely pathogenic": 1,
    "uncertain significance": 2,
    "likely benign": 3, "benign/likely benign": 3,
    "benign": 4, "other": 5, "unknown": 6,
}


def query_clinvar(gene: str, original_aa: str,
                  position: int, new_aa: str) -> dict:
    """
    Search ClinVar for a mutation and return the best matching variant.
    Tries multiple search strategies, prefers simple (non-compound) mutations.
    """
    orig3 = ONE_TO_THREE.get(original_aa, original_aa)
    new3  = ONE_TO_THREE.get(new_aa, new_aa)
    mutation_str = f"{original_aa}{position}{new_aa}"

    print(f"\n   Searching ClinVar for {gene} {mutation_str}...")

    # Strategy 1 — search by protein change (most specific)
    search_strategies = [
        f"{gene}[gene] AND {orig3}{position}{new3}[variant name]",
        f"{gene}[gene] AND p.{orig3}{position}{new3}",
        f"{gene}[gene] AND {original_aa}{position}{new_aa}[variant name]",
        f'"{orig3} {position} {new3}"[All Fields] AND {gene}[gene]',
        f"{gene}[gene]",   # broad fallback
    ]

    all_ids = []
    for term in search_strategies:
        ids = _search_ids(term, retmax=15)
        for i in ids:
            if i not in all_ids:
                all_ids.append(i)
        if len(all_ids) >= 10:
            break

    if not all_ids:
        return {"status": "not_found", "message": f"{gene} {mutation_str} not in ClinVar",
                "gene": gene, "mutation": mutation_str}

    # Fetch all variants
    result = _fetch_variants_by_ids(all_ids[:15])

    # Pick best variant — prefer simple single mutations and high significance
    best_variant  = None
    best_score    = 999

    for vid in all_ids[:15]:
        v = result.get(vid, {})
        if not v:
            continue

        germline   = v.get("germline_classification", {})
        sig        = germline.get("description", "unknown").lower()
        title      = v.get("title", "")
        is_compound = "[" in title and ";" in title

        score = PRIORITY.get(sig, 5) + (3 if is_compound else 0)

        print(f"   [{vid}] sig='{sig}' compound={is_compound} score={score} | {title[:55]}")

        if score < best_score:
            best_score   = score
            best_variant = v

    if not best_variant:
        return {"status": "not_found", "message": f"{gene} {mutation_str} not in ClinVar",
                "gene": gene, "mutation": mutation_str}

    return _parse_variant(best_variant, len(all_ids))


def _parse_variant(variant: dict, total: int) -> dict:
    """Parse a ClinVar variant record into a clean dict."""

    germline      = variant.get("germline_classification", {})
    significance  = germline.get("description", "Unknown")
    review_status = germline.get("review_status", "Unknown")

    # ── Extract and deduplicate diseases ─────────────────────────
    raw_diseases = []

    # Path 1: germline_classification → trait_set
    for trait in germline.get("trait_set", []):
        name = trait.get("trait_name", "")
        if name and name.lower() not in ["not specified", "not provided", "see cases"]:
            raw_diseases.append(name)

    # Path 2: variation_set → trait_set
    for var in variant.get("variation_set", []):
        for trait in var.get("trait_set", []):
            name = trait.get("trait_name", "")
            if name and name.lower() not in ["not specified", "not provided", "see cases"]:
                raw_diseases.append(name)

    # Path 3: supporting submissions
    subs = variant.get("supporting_submissions", {})
    for sub_type in ["germline", "somatic"]:
        for sub in subs.get(sub_type, []):
            for cond in sub.get("condition_list", []):
                name = cond.get("name", "")
                if name and name.lower() not in ["not specified", "not provided", "see cases"]:
                    raw_diseases.append(name)

    # Deduplicate (case-insensitive)
    seen       = set()
    unique_diseases = []
    for d in raw_diseases:
        key = d.lower().strip()
        if key not in seen:
            seen.add(key)
            unique_diseases.append(d)

    # ── Prioritize most relevant disease ─────────────────────────
    # Score each disease — higher score = more relevant
    # We want the most specific, well-known disease name at the top
    def disease_score(name: str) -> int:
        n = name.lower()
        # Penalize vague/broad terms
        if any(x in n for x in ["inborn", "related disorder", "locus", "trait", "susceptibility"]):
            return 10
        # Prefer specific well-known diseases
        if any(x in n for x in ["sickle", "hb ss", "thalassemia", "anemia",
                                  "cancer", "syndrome", "disease", "disorder"]):
            return 0
        return 5

    unique_diseases.sort(key=disease_score)

    # Keep only top 5 most relevant
    diseases = unique_diseases[:5]

    # ── Extract OMIM IDs ──────────────────────────────────────────
    omim_ids = []
    for var in variant.get("variation_set", []):
        for xref in var.get("variation_xrefs", []):
            if xref.get("db_source") == "OMIM":
                base_id = xref.get("db_id", "").split(".")[0]
                if base_id and base_id not in omim_ids:
                    omim_ids.append(base_id)

    variation_id   = str(variant.get("variation_id", "") or variant.get("uid", ""))
    variation_name = variant.get("title", "")

    sig_info = SIGNIFICANCE_MAP.get(
        significance,
        {"level": "unknown", "emoji": "❓", "color": "#95a5a6"}
    )

    print(f"\n   ✅ Final result:")
    print(f"   Significance: {significance}")
    print(f"   Top diseases: {diseases}")
    print(f"   OMIM IDs:     {omim_ids}")

    return {
        "status":         "found",
        "significance":   significance,
        "level":          sig_info["level"],
        "emoji":          sig_info["emoji"],
        "color":          sig_info["color"],
        "diseases":       diseases if diseases else ["Not specified"],
        "review_status":  review_status,
        "variation_name": variation_name,
        "variation_id":   variation_id,
        "total_results":  total,
        "omim_ids":       omim_ids,
        "clinvar_url":    f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/"
                          if variation_id else "",
    }


# ── OMIM lookup using ID directly ─────────────────────────────────────────────

def query_omim_by_id(omim_id: str) -> dict | None:
    """
    Fetch disease info directly from OMIM using the ID we extracted
    from ClinVar's variation_xrefs. No API key needed for this approach
    — we use the NCBI OMIM database via E-utilities.
    """
    try:
        # Search OMIM entry by ID
        resp = _get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "omim", "id": omim_id, "retmode": "json"}
        )
        if not resp or resp.status_code != 200:
            return None

        data  = resp.json()
        entry = data.get("result", {}).get(omim_id, {})

        if not entry:
            return None

        title = entry.get("title", "")
        # OMIM titles are ALL CAPS — convert to Title Case
        title = title.title() if title else ""

        return {
            "omim_id":    omim_id,
            "title":      title,
            "omim_url":   f"https://www.omim.org/entry/{omim_id}",
            "source":     "OMIM via NCBI",
        }

    except Exception as e:
        print(f"   OMIM lookup error: {e}")
        return None


def query_omim(disease_name: str, omim_ids: list = None) -> dict | None:
    """
    Get disease info — tries 3 sources in order:
    1. Wikipedia (best plain English)
    2. NIH MedlinePlus (medical but readable)
    3. NCBI Gene Reviews (technical fallback)
    """

    omim_result = None

    # Get OMIM title via NCBI
    if omim_ids:
        for oid in omim_ids:
            omim_result = query_omim_by_id(oid)
            if omim_result:
                break

    # Try to get description — multiple sources
    description = None
    source_used = None

    # Source 1 — Wikipedia (best for readability)
    search_name = _pick_best_disease_name(disease_name)
    description = _query_wikipedia(search_name)
    if description:
        source_used = "Wikipedia"

    # Source 2 — Try Wikipedia with raw disease name if mapped name failed
    if not description and search_name != disease_name:
        description = _query_wikipedia(disease_name)
        if description:
            source_used = "Wikipedia"

    # Source 3 — NIH MedlinePlus
    if not description:
        description = _query_medlineplus_connect(search_name)
        if description:
            source_used = "MedlinePlus"

    # Source 4 — NCBI Gene database (searches by gene name)
    if not description and omim_ids:
        description = _query_ncbi_gene_description(disease_name)
        if description:
            source_used = "NCBI"

    # Source 5 — Build a basic description from ClinVar data alone
    if not description:
        description = _build_basic_description(disease_name)
        source_used = "ClinVar"

    if omim_result:
        omim_result["description"] = description
        omim_result["source"]      = source_used
        return omim_result

    if description:
        return {
            "title":       search_name,
            "description": description,
            "omim_url":    f"https://www.omim.org/search?search={search_name.replace(' ', '+')}",
            "source":      source_used,
        }

    return None


def _query_ncbi_gene_description(disease_name: str) -> str | None:
    """Search NCBI Gene database for disease description."""
    try:
        # Search for the gene/disease
        search_resp = _get(
            CLINVAR_SEARCH_URL,
            params={
                "db":      "gene",
                "term":    disease_name,
                "retmode": "json",
                "retmax":  1,
            }
        )
        if not search_resp:
            return None

        ids = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None

        # Fetch gene summary
        summary_resp = _get(
            CLINVAR_SUMMARY_URL,
            params={
                "db":      "gene",
                "id":      ids[0],
                "retmode": "json",
            }
        )
        if not summary_resp:
            return None

        result  = summary_resp.json().get("result", {})
        entry   = result.get(ids[0], {})
        summary = entry.get("summary", "")

        if summary and len(summary) > 50:
            return summary[:600]

        return None

    except:
        return None


def _build_basic_description(disease_name: str) -> str:
    """
    Last resort — build a basic description from the disease name alone.
    Better than showing nothing.
    """
    return (
        f"{disease_name} is a genetic condition catalogued in the ClinVar "
        f"and OMIM databases. It is caused by mutations in the associated gene "
        f"that affect protein structure and function. "
        f"Please visit the OMIM and ClinVar links below for full clinical details, "
        f"including symptoms, inheritance pattern, and treatment options."
    )


def _pick_best_disease_name(disease_name: str) -> str:
    """Map ClinVar technical names to common Wikipedia-searchable names."""
    name = disease_name.lower()

    mappings = {
        "hb ss":                      "Sickle cell disease",
        "sickle":                     "Sickle cell disease",
        "beta-thalassemia":           "Beta thalassemia",
        "thalassemia":                "Thalassemia",
        "brca1":                      "BRCA1",
        "brca2":                      "BRCA2",
        "breast":                     "Breast cancer",
        "ovarian":                    "Ovarian cancer",
        "cystic fibrosis":            "Cystic fibrosis",
        "cftr":                       "Cystic fibrosis",
        "parkinson":                  "Parkinson's disease",
        "alzheimer":                  "Alzheimer's disease",
        "huntington":                 "Huntington's disease",
        "hemophilia":                 "Hemophilia",
        "phenylketonuria":            "Phenylketonuria",
        "marfan":                     "Marfan syndrome",
        "lynch":                      "Lynch syndrome",
        "li-fraumeni":                "Li-Fraumeni syndrome",
        "tp53":                       "Li-Fraumeni syndrome",
        "fabry":                      "Fabry disease",
        "gaucher":                    "Gaucher's disease",
        "wilson":                     "Wilson's disease",
        "hemochromatosis":            "Hereditary hemochromatosis",
        "retinoblastoma":             "Retinoblastoma",
        "neurofibromatosis":          "Neurofibromatosis",
        "tuberous sclerosis":         "Tuberous sclerosis",
        "von hippel":                 "Von Hippel-Lindau disease",
        "polycystic kidney":          "Polycystic kidney disease",
        "muscular dystrophy":         "Muscular dystrophy",
        "spinal muscular":            "Spinal muscular atrophy",
        "fragile x":                  "Fragile X syndrome",
        "down syndrome":              "Down syndrome",
        "turner":                     "Turner syndrome",
        "klinefelter":                "Klinefelter syndrome",
        "DiGeorge":                   "DiGeorge syndrome",
        "prader":                     "Prader-Willi syndrome",
        "angelman":                   "Angelman syndrome",
        "rett":                       "Rett syndrome",
        "ataxia":                     "Ataxia",
        "galactosemia":               "Galactosemia",
        "homocystinuria":             "Homocystinuria",
        "maple syrup":                "Maple syrup urine disease",
        "glycogen storage":           "Glycogen storage disease",
        "mucopolysaccharidosis":      "Mucopolysaccharidosis",
        "niemann-pick":               "Niemann-Pick disease",
        "tay-sachs":                  "Tay-Sachs disease",
        "krabbe":                     "Krabbe disease",
        "metachromatic":              "Metachromatic leukodystrophy",
        "adrenoleukodystrophy":       "Adrenoleukodystrophy",
        "hereditary breast":          "Hereditary breast and ovarian cancer",
        "cowden":                     "Cowden syndrome",
        "pten":                       "Cowden syndrome",
        "multiple endocrine":         "Multiple endocrine neoplasia",
        "von willebrand":             "Von Willebrand disease",
        "fanconi":                    "Fanconi anemia",
        "bloom":                      "Bloom syndrome",
        "xeroderma":                  "Xeroderma pigmentosum",
        "ataxia telangiectasia":      "Ataxia-telangiectasia",
        "werner":                     "Werner syndrome",
        "progeria":                   "Progeria",
    }

    for key, mapped in mappings.items():
        if key in name:
            return mapped

    # Generic cleanup — remove database-specific suffixes
    cleaned = disease_name
    for suffix in [" HBB/LCRB", " NM_", ", familial", ", hereditary"]:
        cleaned = cleaned.replace(suffix, "")

    return cleaned.strip()

def _pick_best_disease_name(disease_name: str) -> str:
    """
    Map ClinVar's technical disease names to common searchable names.
    ClinVar uses names like 'Hb SS disease' — Wikipedia knows 'Sickle cell disease'.
    """
    name = disease_name.lower()

    mappings = {
        "hb ss":           "Sickle cell disease",
        "sickle":          "Sickle cell disease",
        "thalassemia":     "Thalassemia",
        "brca":            "Hereditary breast and ovarian cancer",
        "breast":          "Breast cancer",
        "cystic fibrosis": "Cystic fibrosis",
        "cftr":            "Cystic fibrosis",
        "parkinson":       "Parkinson's disease",
        "alzheimer":       "Alzheimer's disease",
        "huntington":      "Huntington's disease",
        "hemophilia":      "Hemophilia",
        "phenylketonuria": "Phenylketonuria",
        "marfan":          "Marfan syndrome",
    }

    for key, mapped in mappings.items():
        if key in name:
            return mapped

    # Return original if no mapping found
    return disease_name


def _query_wikipedia(disease_name: str) -> str | None:
    """
    Fetch plain-English disease description from Wikipedia API.
    Must include User-Agent header — Wikipedia blocks requests without one.
    """
    headers = {
        "User-Agent": "ProteinExplorer/1.0 (educational bioinformatics tool; contact@example.com)"
    }

    try:
        # Try direct title lookup first
        url = "https://en.wikipedia.org/api/rest_v1/page/summary/" + \
              disease_name.replace(" ", "_")

        resp = requests.get(url, headers=headers, timeout=15, verify=False)

        if resp.status_code == 200:
            extract = resp.json().get("extract", "")
            if len(extract) > 50:
                print(f"   ✅ Wikipedia found: '{disease_name}'")
                return extract[:600]

        # Direct lookup failed — try search
        search_resp = requests.get(
            "https://en.wikipedia.org/w/api.php",
            headers=headers,
            params={
                "action":   "query",
                "list":     "search",
                "srsearch": disease_name,
                "format":   "json",
                "srlimit":  1,
            },
            timeout=15,
            verify=False
        )

        if not search_resp or search_resp.status_code != 200:
            return None

        results = search_resp.json().get("query", {}).get("search", [])
        if not results:
            return None

        top_title = results[0]["title"].replace(" ", "_")
        print(f"   Wikipedia search found: '{results[0]['title']}'")

        summary_resp = requests.get(
            f"https://en.wikipedia.org/api/rest_v1/page/summary/{top_title}",
            headers=headers,
            timeout=15,
            verify=False
        )

        if summary_resp and summary_resp.status_code == 200:
            extract = summary_resp.json().get("extract", "")
            if len(extract) > 50:
                return extract[:600]

        return None

    except Exception as e:
        print(f"   Wikipedia error: {e}")
        return None
def _query_medlineplus_connect(disease_name: str) -> str | None:
    """Fallback MedlinePlus Connect API."""
    try:
        resp = _get(
            "https://connect.medlineplus.gov/application",
            params={
                "mainSearchCriteria.v.cs": "2.16.840.1.113883.6.90",
                "knowledgeResponseType":   "application/json",
                "informationRecipient":    "PROV",
                "mainSearchCriteria.v.dn": disease_name,
            }
        )
        if not resp or resp.status_code != 200:
            return None

        entries = resp.json().get("feed", {}).get("entry", [])
        if not entries:
            return None

        summary = entries[0].get("summary", {}).get("_value", "")
        return re.sub(r'<[^>]+>', '', summary).strip()[:600] or None

    except:
        return None


# ── Main entry point ──────────────────────────────────────────────────────────

def get_disease_info(gene: str, original_aa: str,
                     position: int, new_aa: str) -> dict:
    """
    Full pipeline: ClinVar lookup → OMIM lookup → combined result.
    This is what app.py calls.
    """
    result = {
        "clinvar":  None,
        "omim":     None,
        "gene":     gene,
        "mutation": f"{original_aa}{position}{new_aa}",
    }

    # Step 1 — ClinVar
    clinvar = query_clinvar(gene, original_aa, position, new_aa)
    result["clinvar"] = clinvar

    # Step 2 — OMIM (only if ClinVar found something)
    if clinvar.get("status") == "found":
        omim_ids = clinvar.get("omim_ids", [])
        diseases = clinvar.get("diseases", [])
        disease_name = diseases[0] if diseases and diseases[0] != "Not specified" else ""

        print(f"\n   Looking up OMIM (IDs: {omim_ids}, name: '{disease_name}')...")
        omim = query_omim(disease_name, omim_ids)
        result["omim"] = omim

    return result


# ── Test ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("Testing: HBB E6V (Sickle Cell Mutation)")
    print("=" * 60)

    result = get_disease_info("HBB", "E", 6, "V")

    print("\n── ClinVar ──")
    cv = result["clinvar"]
    if cv["status"] == "found":
        print(f"  {cv['emoji']} {cv['significance']}")
        print(f"  Disease:   {cv['diseases']}")
        print(f"  Review:    {cv['review_status']}")
        print(f"  OMIM IDs:  {cv.get('omim_ids', [])}")
        print(f"  URL:       {cv['clinvar_url']}")
    else:
        print(f"  {cv['message']}")

    print("\n── OMIM ──")
    om = result["omim"]
    if om:
        print(f"  Title:       {om.get('title', 'N/A')}")
        print(f"  Description: {str(om.get('description', 'N/A'))[:200]}...")
        print(f"  URL:         {om.get('omim_url', 'N/A')}")
    else:
        print("  Not found")