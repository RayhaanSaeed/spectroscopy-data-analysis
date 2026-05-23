import os
import csv
import re
import math
import json
from datetime import datetime
from collections import defaultdict
from urllib.parse import urlencode
from urllib.request import urlopen, Request

import chardet

# -----------------------------
# CONFIG (portable)
# -----------------------------
PROJECT_FOLDER_NAME = "spectroscopt db project"
CONVERTED_SUBFOLDER = os.path.join("Project csv data", "converted_opju")
MASTER_NAME = "master_spectra_export.csv"

MAX_PER_PROTEIN = 2              # Absorption + Emission (ideal)
MIN_POINTS_NUMERIC = 40          # numeric points threshold (works for FPbase & most papers)
MIN_NONZERO_POINTS = 20          # nonzero numeric points threshold

# When FPbase provides "ex", treat it as Absorption so we get 2 rows (Abs + Em).
FPBASE_EX_COUNTS_AS = "Absorption"

# -----------------------------
# Helpers: paths/encoding/CSV parsing
# -----------------------------
def find_project_root():
    """Find the project root under OneDrive - University of Leeds for any username."""
    home = os.path.expanduser("~")
    candidates = [
        os.path.join(home, "OneDrive - University of Leeds", PROJECT_FOLDER_NAME),
        os.path.join(home, "OneDrive - University of Leeds", PROJECT_FOLDER_NAME.lower()),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c

    onedrive_root = os.path.join(home, "OneDrive - University of Leeds")
    if os.path.exists(onedrive_root):
        for root, dirs, _ in os.walk(onedrive_root):
            for d in dirs:
                if d.lower() == PROJECT_FOLDER_NAME.lower():
                    return os.path.join(root, d)
    return None


def detect_encoding(path):
    with open(path, "rb") as f:
        raw = f.read(4096)
    r = chardet.detect(raw)
    return r.get("encoding") or "utf-8"


def sniff_delimiter(path, encoding):
    """
    Robust delimiter detection:
    - Try csv.Sniffer on a small sample.
    - Fallback to tab if sample contains many tabs, else comma, else semicolon.
    Never crashes.
    """
    try:
        with open(path, "r", encoding=encoding, errors="replace") as f:
            sample = f.read(4096)
        # quick heuristics first (Sniffer can fail on messy samples)
        if sample.count("\t") > sample.count(",") and sample.count("\t") > sample.count(";"):
            return "\t"
        if sample.count(",") >= sample.count(";"):
            return ","
        return ";"
    except Exception:
        return ","


def safe_float(v):
    try:
        if v is None:
            return math.nan
        s = str(v).strip()
        if s in ("", "--", "—", "N/A", "na", "NaN", "nan"):
            return math.nan
        return float(s)
    except Exception:
        return math.nan


def is_usable_series(y):
    """Reject empty / almost-empty / flat series."""
    yy = [v for v in y if not math.isnan(v)]
    if len(yy) < MIN_POINTS_NUMERIC:
        return False
    nonzero = sum(1 for v in yy if abs(v) > 0)
    if nonzero < MIN_NONZERO_POINTS:
        return False
    if (max(yy) - min(yy)) < 1e-12:
        return False
    return True


def normalize_0_1(y):
    yy = [v for v in y if not math.isnan(v)]
    if not yy:
        return None
    lo, hi = min(yy), max(yy)
    if abs(hi - lo) < 1e-12:
        return None
    out = []
    for v in y:
        if math.isnan(v):
            out.append(math.nan)
        else:
            out.append((v - lo) / (hi - lo))
    return out


def peak_and_fwhm(x, y):
    pairs = [(xx, yy) for xx, yy in zip(x, y) if not math.isnan(xx) and not math.isnan(yy)]
    if len(pairs) < 10:
        return None, None
    xs, ys = zip(*pairs)
    maxy = max(ys)
    if maxy <= 0:
        return None, None
    i0 = ys.index(maxy)
    peak = xs[i0]
    half = maxy / 2.0
    above = [xx for xx, yy in zip(xs, ys) if yy >= half]
    if len(above) < 2:
        return peak, None
    return peak, (above[-1] - above[0])


# -----------------------------
# Protein/spectrum classification + filtering
# -----------------------------
NON_PROTEIN_TOKENS = {
    "time", "t", "au", "a.u.", "a.u", "abs", "absorbance", "em", "ex",
    "nm", "wavelength", "intensity", "counts", "signal", "baseline", "fit",
    "residual", "resid", "error", "tau", "lifetime", "decay"
}


def classify(col_name):
    s = (col_name or "").lower()
    if "overlap" in s or "fret" in s:
        return "Ignore"
    if "abs" in s or "absorption" in s:
        return "Absorption"
    # emission: keep first to avoid "em" inside other words being weird
    if " emission" in s or s.endswith(" em") or "emis" in s or "fluor" in s:
        return "Emission"
    if " ex" in s or "excitation" in s or s.endswith(" ex"):
        return "Excitation"
    return "Unknown"


def looks_like_fpbase_header(header_cells):
    """
    FPbase downloads typically look like:
      wavelength    EBFP2 2p   EBFP2 em   EBFP2 ex
    """
    if not header_cells or len(header_cells) < 2:
        return False
    h0 = (header_cells[0] or "").strip().lower()
    if h0 not in {"wavelength", "wl", "lambda"}:
        return False
    joined = " ".join((c or "").lower() for c in header_cells[1:])
    # very FPbase-ish
    return (" em" in joined) or (" ex" in joined) or (" 2p" in joined)


def clean_protein_token(token):
    token = (token or "").strip()
    token = re.sub(r"\s+", " ", token)
    token = token.replace("_", " ").strip()
    return token


def extract_protein_from_yname(yname):
    """
    Stronger protein-name extraction:
    - Take leading part before known suffixes like ' em', ' ex', ' abs', ' absorption'
    - Remove trailing tokens like '2p'
    - Return None if it becomes junk.
    """
    s = clean_protein_token(yname)
    low = s.lower()

    # strip common suffixes
    for suf in [" absorption", " absorbance", " abs", " emission", " em", " excitation", " ex"]:
        if low.endswith(suf):
            s = s[: -len(suf)].strip()
            low = s.lower()

    # remove "2p" token at end
    s = re.sub(r"\b2p\b$", "", s, flags=re.IGNORECASE).strip()

    # if still multi tokens, keep first token if it looks like a FP/protein symbol, else keep full
    s = s.strip()

    return s if protein_name_is_valid(s) else None


def protein_name_is_valid(name):
    """
    Accepts real protein-like names, rejects:
      - single letters (a, b, c)
      - numeric-only
      - known non-protein tokens
      - very short garbage
    """
    if not name:
        return False
    n = clean_protein_token(name)
    low = n.lower()

    if low in NON_PROTEIN_TOKENS:
        return False
    if len(n) <= 1:
        return False
    if re.fullmatch(r"[a-zA-Z]", n):
        return False
    if re.fullmatch(r"\d+", n):
        return False
    # must contain at least one letter
    if not re.search(r"[A-Za-z]", n):
        return False
    # reject pure units-ish
    if low in {"ms", "s", "sec", "min", "hour", "k", "c", "m", "mm", "um"}:
        return False
    return True


# Manual hints for common paper proteins that aren't FPbase proteins
def protein_guess_from_paper(col_name):
    s = (col_name or "").lower()
    if "lhcii" in s or "lhc2" in s or "light-harvesting complex" in s or re.search(r"\blhc\b", s):
        return "LHCII"
    if "texas red" in s or "texas" in s:
        return "Texas Red"
    if "gfp" in s:
        return "GFP"
    return None


# -----------------------------
# Paper PDF metadata extraction (basic but works)
# -----------------------------
def find_pdf_in_folder(folder):
    try:
        for f in os.listdir(folder):
            if f.lower().endswith(".pdf"):
                return os.path.join(folder, f)
    except Exception:
        pass
    return None


def extract_pdf_basic_meta(pdf_path):
    """
    Lightweight: DOI, Year, Title, Authors (+ try journal/volume/pages, pH/temp/conc/exc if present).
    If PyPDF2 isn't installed, returns N/A.
    """
    meta = {
        "SourceType": "Paper",
        "SourceName": os.path.basename(pdf_path) if pdf_path else "N/A",
        "DOI": "N/A",
        "Year": "N/A",
        "Title": "N/A",
        "Authors": "N/A",
        "Journal": "N/A",
        "Volume": "N/A",
        "Pages": "N/A",
        "pH": "N/A",
        "Temperature_K": "N/A",
        "Concentration_M": "N/A",
        "Instrument": "N/A",
        "ExcitationWavelength_nm": "N/A",
        "Category": "N/A",
        "ChromophoreType": "N/A",
        "ProteinClass": "N/A",
        "Credit": "N/A",
    }
    try:
        import PyPDF2
    except Exception:
        return meta
    if not pdf_path or not os.path.exists(pdf_path):
        return meta

    try:
        reader = PyPDF2.PdfReader(pdf_path)
        text = ""
        for p in reader.pages[:4]:
            text += "\n" + (p.extract_text() or "")
        low = text.lower()

        m = re.search(r"10\.\d{4,9}/\S+", text)
        if m:
            meta["DOI"] = m.group(0).rstrip(").,;")
        y = re.search(r"(20\d{2}|19\d{2})", text)
        if y:
            meta["Year"] = y.group(1)

        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if lines:
            meta["Title"] = lines[0][:250]
        for ln in lines[1:12]:
            if "," in ln and len(ln.split()) <= 35:
                meta["Authors"] = ln[:250]
                break

        # crude journal/vol/pages patterns
        jp = re.search(r"([A-Za-z][A-Za-z \-&]+)\s+(\d+)\s*,\s*([0-9]+[-–][0-9]+)\s*\((20\d{2}|19\d{2})\)", text)
        if jp:
            meta["Journal"] = jp.group(1).strip()[:250]
            meta["Volume"]  = jp.group(2).strip()
            meta["Pages"]   = jp.group(3).replace(" ", "")

        ph = re.search(r"pH\s*=?\s*(\d+(\.\d+)?)", text, re.IGNORECASE)
        if ph:
            meta["pH"] = ph.group(1)

        t = re.search(r"(\d+(\.\d+)?)\s*°\s*C", text)
        if t:
            meta["Temperature_K"] = f"{float(t.group(1)) + 273.15:.2f}"

        conc = re.search(r"(\d+(\.\d+)?)\s*(mM|µM|uM|μM|nM|M)\b", text)
        if conc:
            val = float(conc.group(1))
            unit = conc.group(3)
            factor = {"M":1.0,"mM":1e-3,"µM":1e-6,"uM":1e-6,"μM":1e-6,"nM":1e-9}.get(unit, 1.0)
            meta["Concentration_M"] = f"{val*factor:.3e}"

        exc = re.search(r"excitation\s+(?:at|wavelength)\s*(\d+)\s*nm", text, re.IGNORECASE)
        if exc:
            meta["ExcitationWavelength_nm"] = exc.group(1)

        if "fluorometer" in low:
            meta["Instrument"] = "Fluorometer"
        elif "spectrophotometer" in low:
            meta["Instrument"] = "Spectrophotometer"

        # quick category hints
        if "light-harvesting complex" in low or "lhcii" in low or "photosystem" in low:
            meta["Category"] = "Photosynthetic protein complex"
            meta["ProteinClass"] = "Photosynthetic complex"
        if "texas red" in low:
            meta["ChromophoreType"] = "Texas Red"
        if meta["Authors"] != "N/A":
            meta["Credit"] = meta["Authors"]

    except Exception:
        pass

    return meta


# -----------------------------
# FPbase API enrichment (metadata source)
# -----------------------------
FPBASE_PROTEINS_ENDPOINT = "https://www.fpbase.org/api/proteins/"

def _http_get_json(url: str, timeout_s: int = 25) -> dict:
    req = Request(
        url,
        headers={"User-Agent": "spectroscopy-db-importer/1.0", "Accept": "application/json"},
    )
    with urlopen(req, timeout=timeout_s) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    return json.loads(raw)

def fpbase_fetch_protein_record(protein_name: str, timeout_s: int = 25):
    name = (protein_name or "").strip()
    if not name:
        return None

    def build_url(params: dict) -> str:
        return FPBASE_PROTEINS_ENDPOINT + "?" + urlencode(params)

    # exact match then contains
    for params in [{"format":"json","name__iexact":name}, {"format":"json","name__icontains":name}]:
        try:
            data = _http_get_json(build_url(params), timeout_s=timeout_s)
            results = data.get("results", []) if isinstance(data, dict) else []
            if results:
                # prefer exact name
                nl = name.lower()
                for r in results:
                    if str(r.get("name","")).lower() == nl:
                        return r
                return results[0]
        except Exception:
            continue
    return None

def fpbase_extract_meta(fpbase_record: dict, protein_name: str):
    """
    Returns a meta dict aligned to your DB fields.
    (Fail-safe: returns N/A if missing)
    """
    meta = {
        "SourceType": "FPbase",
        "SourceName": "FPbase",
        "DOI": "N/A",
        "Year": "N/A",
        "Title": "N/A",
        "Authors": "N/A",
        "Journal": "N/A",
        "Volume": "N/A",
        "Pages": "N/A",
        "pH": "N/A",
        "Temperature_K": "N/A",
        "Concentration_M": "N/A",
        "Instrument": "N/A",
        "ExcitationWavelength_nm": "N/A",
        "Category": "Fluorescent protein",
        "ChromophoreType": "N/A",
        "ProteinClass": "FP",
        "Credit": "FPbase",
        "UniProtID": "N/A",
        "PDB_ID": "N/A",
        "StructureLink": "N/A",
        "FPbase_URL": "N/A",
    }

    if not fpbase_record:
        return meta

    slug = fpbase_record.get("slug") or ""
    if slug:
        meta["FPbase_URL"] = f"https://www.fpbase.org/protein/{slug}/"

    # FPbase page has UniProtKB link; API may include it (not guaranteed). Try common keys:
    for k in ["uniprot", "uniprot_id", "uniprotkb", "uniprot_accession", "uniprot_accession_id"]:
        v = fpbase_record.get(k)
        if v and isinstance(v, str):
            meta["UniProtID"] = v.strip()
            break

    # Primary reference (if present) may have DOI, citation etc.
    # API schemas can differ; keep best-effort.
    prim = fpbase_record.get("primary_reference") or fpbase_record.get("primaryReference") or None
    if isinstance(prim, dict):
        meta["Title"] = prim.get("title") or meta["Title"]
        meta["Year"] = str(prim.get("year") or meta["Year"])
        meta["DOI"] = prim.get("doi") or meta["DOI"]
        meta["Journal"] = prim.get("journal") or meta["Journal"]
        meta["Volume"] = str(prim.get("volume") or meta["Volume"])
        meta["Pages"] = prim.get("pages") or meta["Pages"]
        auth = prim.get("authors")
        if auth:
            meta["Authors"] = auth if isinstance(auth, str) else ", ".join(map(str, auth))
        meta["Credit"] = "FPbase"

    # Chromophore/cofactor sometimes present
    for k in ["cofactor", "chromophore", "chromophore_type"]:
        v = fpbase_record.get(k)
        if v:
            meta["ChromophoreType"] = str(v)[:200]
            break

    return meta


# -----------------------------
# UniProt + PDB enrichment (best effort)
# -----------------------------
_UNIPROT_CACHE = {}
_PDB_CACHE = {}

def uniprot_lookup_accession(protein_name: str, timeout_s: int = 25):
    """
    Best-effort UniProt accession lookup by protein name.
    Uses UniProt REST search. Fail-safe to None.
    """
    name = (protein_name or "").strip()
    if not name:
        return None
    if name in _UNIPROT_CACHE:
        return _UNIPROT_CACHE[name]

    try:
        # UniProt REST: query by protein name in recommendedName OR gene/protein name.
        # Request a small JSON response.
        base = "https://rest.uniprot.org/uniprotkb/search"
        query = f'(protein_name:"{name}" OR gene:"{name}" OR id:{name}) AND reviewed:true'
        url = base + "?" + urlencode({"query": query, "format": "json", "size": 1})
        data = _http_get_json(url, timeout_s=timeout_s)
        results = data.get("results", [])
        if results:
            acc = results[0].get("primaryAccession")
            _UNIPROT_CACHE[name] = acc
            return acc
    except Exception:
        pass

    _UNIPROT_CACHE[name] = None
    return None


def rcsb_best_pdb_for_uniprot(uniprot_acc: str, timeout_s: int = 25):
    """
    Best-effort: find a PDB ID for a UniProt accession via RCSB search API.
    """
    if not uniprot_acc:
        return None
    if uniprot_acc in _PDB_CACHE:
        return _PDB_CACHE[uniprot_acc]

    try:
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "value": uniprot_acc
                }
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 1}}
        }
        req = Request(url, data=json.dumps(payload).encode("utf-8"),
                      headers={"Content-Type": "application/json", "Accept": "application/json"})
        with urlopen(req, timeout=timeout_s) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
        hits = data.get("result_set", [])
        if hits:
            pdb = hits[0].get("identifier")
            _PDB_CACHE[uniprot_acc] = pdb
            return pdb
    except Exception:
        pass

    _PDB_CACHE[uniprot_acc] = None
    return None


def enrich_uniprot_pdb(meta: dict, protein_name: str):
    """
    Fill UniProtID/PDB_ID/StructureLink where possible.
    Only attempts for FPbase proteins (or protein-like names).
    """
    if meta.get("SourceType") != "FPbase":
        # Don't force this for paper complexes like LHCII
        return meta

    unip = meta.get("UniProtID")
    if not unip or unip == "N/A":
        unip = uniprot_lookup_accession(protein_name)
        if unip:
            meta["UniProtID"] = unip

    pdb = None
    if unip and unip != "N/A":
        pdb = rcsb_best_pdb_for_uniprot(unip)
        if pdb:
            meta["PDB_ID"] = pdb
            meta["StructureLink"] = f"https://www.rcsb.org/structure/{pdb}"

    return meta


# -----------------------------
# Core: scan CSVs and select
# -----------------------------
def scan_csv(csv_path, paper_meta):
    """
    Returns candidate spectra dicts.
    Applies strict filtering:
      - Reject non-protein series by name
      - Require usable y-series
      - Route metadata correctly (paper vs FPbase)
    """
    enc = detect_encoding(csv_path)
    delim = sniff_delimiter(csv_path, enc)

    with open(csv_path, "r", encoding=enc, errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        rows = list(reader)

    if not rows or len(rows[0]) < 2:
        return []

    header = rows[0]
    data = rows[1:]

    # decide if this CSV is FPbase-like
    is_fpbase = looks_like_fpbase_header(header) or ("fpbase" in os.path.basename(csv_path).lower())

    x_name = header[0]
    y_names = header[1:]

    # x column
    x = []
    for r in data:
        if len(r) >= 1:
            x.append(safe_float(r[0]))
        else:
            x.append(math.nan)

    out = []

    for j, yname in enumerate(y_names, start=1):
        raw_type = classify(yname)
        if raw_type == "Ignore":
            continue

        # y values
        y = []
        for r in data:
            if len(r) > j:
                y.append(safe_float(r[j]))
            else:
                y.append(math.nan)

        if not is_usable_series(y):
            continue

        peak, fwhm = peak_and_fwhm(x, y)
        if peak is None:
            continue

        # protein extraction
        prot = None
        if is_fpbase:
            prot = extract_protein_from_yname(yname)
        else:
            prot = protein_guess_from_paper(yname) or extract_protein_from_yname(yname)

        # If still unknown -> reject (prevents A/B/C/time/AU junk)
        if not protein_name_is_valid(prot):
            continue

        # spectrum type mapping
        spec_type = raw_type
        if is_fpbase and spec_type == "Excitation":
            spec_type = FPBASE_EX_COUNTS_AS  # treat FPbase ex as Absorption

        # dynamic range
        yy = [v for v in y if not math.isnan(v)]
        dyn = max(yy) - min(yy)

        # metadata: paper vs fpbase
        meta = None
        if is_fpbase:
            rec = fpbase_fetch_protein_record(prot)
            meta = fpbase_extract_meta(rec, prot)
            meta = enrich_uniprot_pdb(meta, prot)
        else:
            meta = dict(paper_meta)  # copy

        out.append({
            "protein": prot,
            "spectrum_type": spec_type,
            "csv_path": csv_path,
            "csv": os.path.basename(csv_path),
            "y_name": yname,
            "x_name": x_name,
            "x": x,
            "y": y,
            "peak": peak,
            "fwhm": fwhm,
            "dyn_range": dyn,
            "meta": meta,
            "is_fpbase": is_fpbase,
        })

    return out


def select_best(candidates):
    """
    For each protein:
      keep best Absorption + best Emission (by dyn_range).
      cap at MAX_PER_PROTEIN.
    """
    by_prot = defaultdict(list)
    for c in candidates:
        by_prot[c["protein"]].append(c)

    selected = []
    for prot, specs in by_prot.items():
        best_abs = None
        best_em = None

        for s in specs:
            if s["spectrum_type"] == "Absorption":
                if best_abs is None or s["dyn_range"] > best_abs["dyn_range"]:
                    best_abs = s
            elif s["spectrum_type"] == "Emission":
                if best_em is None or s["dyn_range"] > best_em["dyn_range"]:
                    best_em = s

        picks = []
        if best_abs:
            picks.append(best_abs)
        if best_em:
            picks.append(best_em)

        if not picks:
            # fallback: pick the single strongest
            strongest = max(specs, key=lambda z: z["dyn_range"])
            picks = [strongest]

        picks = picks[:MAX_PER_PROTEIN]
        selected.extend(picks)

    return selected


def write_master_csv(out_path, selected):
    """
    Writes a MASTER CSV with ALL the fields you agreed for DB import,
    plus raw+normalized arrays for Origin plotting.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # compute stokes per protein (using selected abs/em peaks)
    by_prot_type = defaultdict(dict)
    for s in selected:
        by_prot_type[s["protein"]][s["spectrum_type"]] = s

    stokes = {}
    for prot, m in by_prot_type.items():
        if "Absorption" in m and "Emission" in m:
            stokes[prot] = m["Emission"]["peak"] - m["Absorption"]["peak"]

    # Full DB-aligned fields + origin-oriented arrays
    fields = [
        "MoleculeName","ChromophoreType","Category","UniProtID","PDB_ID","StructureLink","ProteinClass",
        "DOI","Authors","Year","Title","Journal","Volume","Pages","Credit",
        "ImportSource","ImportDate","SpectrumType","FileFormat",
        "X_ColumnName","Y_ColumnName","X_Units","Y_Units",
        "X_Values","Y_Values","DataFilePath","Solvent","pH","Temperature_K","Concentration_M",
        "Instrument","ExcitationWavelength_nm","Notes","PeakWavelength_nm","FWHM_nm","StokesShift_nm",

        # extra provenance fields (do NOT break Access; these are for master CSV usefulness)
        "SourceType","SourceName","FPbase_URL",

        # origin plotting convenience
        "X_raw","Y_raw","Y_norm"
    ]

    def pack(arr):
        return ",".join("" if math.isnan(v) else str(v) for v in arr)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()

        for s in selected:
            yn = normalize_0_1(s["y"])
            if yn is None:
                continue

            meta = s["meta"] or {}
            source_type = meta.get("SourceType", "N/A")
            source_name = meta.get("SourceName", "N/A")

            # molecule naming: clear & stable
            # Keep it simple: Protein – SpectrumType (SourceCSV / YColumnName)
            mol_name = f"{s['protein']} – {s['spectrum_type']} ({s['csv']} / {s['y_name']})"

            row = {k: "N/A" for k in fields}

            # DB-aligned core
            row["MoleculeName"] = mol_name
            row["SpectrumType"] = s["spectrum_type"]
            row["FileFormat"] = "csv"
            row["ImportSource"] = os.path.basename(s["csv"])
            row["ImportDate"] = datetime.now().strftime("%Y-%m-%d")

            row["X_ColumnName"] = s["x_name"]
            row["Y_ColumnName"] = s["y_name"]
            row["X_Units"] = "nm"
            row["Y_Units"] = "a.u."

            # Access-friendly fields: keep existing X_Values/Y_Values pattern
            row["X_Values"] = pack(s["x"])
            row["Y_Values"] = pack(s["y"])
            row["DataFilePath"] = s["csv_path"]

            row["PeakWavelength_nm"] = s["peak"]
            row["FWHM_nm"] = s["fwhm"] if s["fwhm"] is not None else ""
            row["StokesShift_nm"] = stokes.get(s["protein"], "")

            # Meta (correct routing already handled)
            row["DOI"] = meta.get("DOI", "N/A")
            row["Authors"] = meta.get("Authors", "N/A")
            row["Year"] = meta.get("Year", "N/A")
            row["Title"] = meta.get("Title", "N/A")
            row["Journal"] = meta.get("Journal", "N/A")
            row["Volume"] = meta.get("Volume", "N/A")
            row["Pages"] = meta.get("Pages", "N/A")
            row["Credit"] = meta.get("Credit", "N/A")
            row["Category"] = meta.get("Category", "N/A")
            row["ChromophoreType"] = meta.get("ChromophoreType", "N/A")
            row["ProteinClass"] = meta.get("ProteinClass", "N/A")

            row["pH"] = meta.get("pH", "N/A")
            row["Temperature_K"] = meta.get("Temperature_K", "N/A")
            row["Concentration_M"] = meta.get("Concentration_M", "N/A")
            row["Instrument"] = meta.get("Instrument", "N/A")
            row["ExcitationWavelength_nm"] = meta.get("ExcitationWavelength_nm", "N/A")

            row["UniProtID"] = meta.get("UniProtID", "N/A")
            row["PDB_ID"] = meta.get("PDB_ID", "N/A")
            row["StructureLink"] = meta.get("StructureLink", "N/A")

            # Provenance extras
            row["SourceType"] = source_type
            row["SourceName"] = source_name
            row["FPbase_URL"] = meta.get("FPbase_URL", "N/A")

            # Notes
            row["Notes"] = f"Selected best spectrum for {s['protein']} / {s['spectrum_type']}"

            # Origin arrays
            row["X_raw"] = pack(s["x"])
            row["Y_raw"] = pack(s["y"])
            row["Y_norm"] = pack(yn)

            w.writerow(row)


# -----------------------------
# MAIN
# -----------------------------
def main():
    project_root = find_project_root()
    if not project_root:
        print("❌ Could not find project root automatically.")
        return

    data_folder = os.path.join(project_root, CONVERTED_SUBFOLDER)
    master_path = os.path.join(project_root, "Project csv data", MASTER_NAME)

    print("=== Building master spectra CSV (selective, portable, correct metadata) ===")
    print("PROJECT_ROOT:", project_root)
    print("DATA_FOLDER:", data_folder)
    print("MASTER_OUT:", master_path)

    if not os.path.exists(data_folder):
        print("❌ DATA_FOLDER does not exist:", data_folder)
        return

    pdf_path = find_pdf_in_folder(data_folder)
    paper_meta = extract_pdf_basic_meta(pdf_path)
    print("PDF:", os.path.basename(pdf_path) if pdf_path else "None")
    print("PDF_META DOI:", paper_meta.get("DOI", "N/A"))

    csv_files = [os.path.join(data_folder, f) for f in os.listdir(data_folder) if f.lower().endswith(".csv")]
    print(f"Found {len(csv_files)} CSV files.")

    all_candidates = []
    for c in csv_files:
        try:
            all_candidates.extend(scan_csv(c, paper_meta))
        except Exception as e:
            print(f"  ⚠ Failed scanning {os.path.basename(c)}: {e}")

    print("Total usable spectra candidates found:", len(all_candidates))

    selected = select_best(all_candidates)
    print("Selected spectra:", len(selected))

    if not selected:
        print("⚠ Selected 0 spectra. Master CSV would be empty.")
        print("   This means columns were filtered out or protein names were invalid.")
        return

    write_master_csv(master_path, selected)
    print("\n✅ Wrote master CSV:", master_path)
    print("It should now contain: 2 rows per protein (Abs + Em) when available,")
    print("with correct per-source metadata and UniProt/PDB best-effort for FPbase proteins.")


if __name__ == "__main__":
    main()
