# ============================================================
# Master CSV  →  Origin Exporter (NO REPEAT + 2 GRAPHS + OVERLAYS)
# - Reads: Project csv data/master_spectra_export.csv
# - Writes per-protein workbook from TSV import (robust)
# - Creates two graphs per protein: Absorption + Emission
# - Creates batch overlay workbooks/graphs for easy comparison
# - Uses a progress JSON so batches never repeat
# ============================================================

import os
import csv
import json
import math
import time
from collections import defaultdict

import win32com.client as win32


# =========================
# CONFIG
# =========================
PROJECT_FOLDER_NAME = "spectroscopt db project"
MASTER_REL_PATH = os.path.join("Project csv data", "master_spectra_export.csv")

# Auto-progress file so re-running continues where it left off
PROGRESS_REL_PATH = os.path.join("Project csv data", "_export_progress_mastercsv.json")

EXPORT_NAME_BASE = "SpectraExport_mastercsv"

BATCH_SIZE = 15            # how many proteins per run
PLOT_NORMALISED = True     # True = plot norm columns, False = raw columns

# Overlay grid (nm). Keep simple.
OVERLAY_GRID_STEP_NM = 1.0


# =========================
# PATH HELPERS
# =========================
def find_project_root():
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


def lt_escape_path(p):
    return p.replace("\\", "\\\\")


def safe_page_name(name, max_len=25):
    allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"
    cleaned = "".join(ch if ch in allowed else "_" for ch in (name or "Protein"))
    cleaned = cleaned.strip("_") or "Protein"
    return cleaned[:max_len]


# =========================
# ROBUST CSV READING
# =========================
def sniff_delimiter(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)

    if sample.count("\t") > sample.count(",") and sample.count("\t") > sample.count(";"):
        return "\t"

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";"])
        return dialect.delimiter
    except Exception:
        # fallback
        if "\t" in sample:
            return "\t"
        if "," in sample:
            return ","
        if ";" in sample:
            return ";"
        return ","


def read_master_rows(master_csv):
    delim = sniff_delimiter(master_csv)
    rows = []
    with open(master_csv, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        for r in reader:
            if r:
                rows.append(r)
    return rows


# =========================
# SERIES PARSER + NORMALISATION
# =========================
def parse_series(cell):
    if cell is None:
        return None
    if not isinstance(cell, str):
        cell = str(cell)
    s = cell.strip()
    if not s:
        return None

    out = []
    for part in s.split(","):
        part = part.strip()
        if not part or part in ("--", "—", "N/A", "na", "NaN", "nan"):
            continue
        try:
            v = float(part)
            if math.isfinite(v):
                out.append(v)
        except Exception:
            continue
    return out if len(out) >= 5 else None


def normalize_0_1(y):
    yy = [v for v in y if math.isfinite(v)]
    if len(yy) < 5:
        return None
    lo, hi = min(yy), max(yy)
    if abs(hi - lo) < 1e-15:
        return None
    return [(v - lo) / (hi - lo) if math.isfinite(v) else float("nan") for v in y]


def classify_type(spectrum_type, y_colname=""):
    st = (spectrum_type or "").lower()
    yc = (y_colname or "").lower()

    # IMPORTANT: treat excitation as Absorption bucket
    if "abs" in st or "absorption" in st or "ex" in st or "excitation" in st:
        return "Absorption"
    if "em" in st or "emis" in st or "fluor" in st:
        return "Emission"

    if "abs" in yc or "ex" in yc:
        return "Absorption"
    if "em" in yc:
        return "Emission"
    return "Other"


def load_spectra_from_master(master_csv):
    """
    Returns: dict protein -> {'Absorption': spec, 'Emission': spec}
    spec: {'x', 'y_raw', 'y_norm'}
    Keeps best (largest dynamic range) per type.
    """
    rows = read_master_rows(master_csv)
    proteins = defaultdict(dict)

    for r in rows:
        protein = (r.get("Protein") or r.get("MoleculeName") or "").strip()
        if not protein:
            continue

        stype = classify_type(r.get("SpectrumType", ""), r.get("YColumnName", ""))
        if stype not in ("Absorption", "Emission"):
            continue

        x = parse_series(r.get("X_raw") or r.get("X_Values") or "")
        y_raw = parse_series(r.get("Y_raw") or r.get("Y_Values") or "")
        if not x or not y_raw:
            continue

        n = min(len(x), len(y_raw))
        if n < 5:
            continue
        x = x[:n]
        y_raw = y_raw[:n]

        y_norm = parse_series(r.get("Y_norm") or "")
        if y_norm:
            y_norm = y_norm[:n]
        else:
            y_norm = normalize_0_1(y_raw)

        if not y_norm or len(y_norm) < 5:
            continue

        dyn = max(y_raw) - min(y_raw)

        existing = proteins[protein].get(stype)
        if existing is None or dyn > existing.get("_dyn", -1):
            proteins[protein][stype] = {"x": x, "y_raw": y_raw, "y_norm": y_norm, "_dyn": dyn}

    # cleanup
    for p in list(proteins.keys()):
        for t in list(proteins[p].keys()):
            proteins[p][t].pop("_dyn", None)
        if not proteins[p]:
            proteins.pop(p, None)

    return proteins


# =========================
# PROGRESS (NO REPEAT)
# =========================
def load_progress(progress_path):
    if not os.path.exists(progress_path):
        return {"exported_proteins": []}
    try:
        with open(progress_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {"exported_proteins": []}


def save_progress(progress_path, progress_obj):
    os.makedirs(os.path.dirname(progress_path), exist_ok=True)
    with open(progress_path, "w", encoding="utf-8") as f:
        json.dump(progress_obj, f, indent=2)


# =========================
# ORIGIN COM HELPERS
# =========================
def get_origin_instance(visible=True):
    # Create a NEW instance (more reliable than attach when you're repeatedly exporting)
    origin = win32.DispatchEx("Origin.ApplicationSI")
    origin.Visible = visible
    return origin


def write_temp_tsv(temp_path, x, abs_raw, abs_norm, em_raw, em_norm):
    os.makedirs(os.path.dirname(temp_path), exist_ok=True)

    def fmt(v):
        try:
            if v is None:
                return ""
            if isinstance(v, float) and (math.isnan(v) or not math.isfinite(v)):
                return ""
        except Exception:
            pass
        return str(v)

    with open(temp_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Wavelength_nm", "Absorption_raw", "Absorption_norm", "Emission_raw", "Emission_norm"])
        for i in range(len(x)):
            w.writerow([fmt(x[i]), fmt(abs_raw[i]), fmt(abs_norm[i]), fmt(em_raw[i]), fmt(em_norm[i])])


def create_workbook_and_import_tsv(origin, protein_name, tsv_path):
    short = safe_page_name(protein_name)
    longn = (protein_name or short).replace('"', "'")

    origin.Execute('newbook name:=TmpBook sheet:=1 option:=lsname;')
    origin.Execute(f'page.name$ = "{short}";')
    origin.Execute(f'page.longname$ = "{longn}";')

    lt_path = lt_escape_path(tsv_path)
    origin.Execute(f'impASC fname:="{lt_path}" options:="TNAME:0,DELIM:9,ONAME:1";')

    # enforce column names
    origin.Execute('wks.col1.lname$="Wavelength"; wks.col1.unit$="nm";')
    origin.Execute('wks.col2.lname$="Absorption_raw"; wks.col2.unit$="a.u.";')
    origin.Execute('wks.col3.lname$="Absorption_norm"; wks.col3.unit$="norm";')
    origin.Execute('wks.col4.lname$="Emission_raw"; wks.col4.unit$="a.u.";')
    origin.Execute('wks.col5.lname$="Emission_norm"; wks.col5.unit$="norm";')

    return short


def create_single_graph(origin, book_name, protein_name, which="Absorption", plot_norm=True):
    """
    which: "Absorption" or "Emission"
    """
    origin.Execute(f'page.name$ = "{book_name}";')

    if which == "Absorption":
        ycol = 3 if plot_norm else 2
        ylab = "Absorption (norm)" if plot_norm else "Absorption (a.u.)"
        title = f"{protein_name} – Absorption"
    else:
        ycol = 5 if plot_norm else 4
        ylab = "Emission (norm)" if plot_norm else "Emission (a.u.)"
        title = f"{protein_name} – Emission"

    # If the entire column is blank/nan, skip to avoid blank graphs.
    # (Origin won't easily tell us; we just allow plotting and you’ll see if it’s empty.
    # But we prevent the earlier bug: emission accidentally plotted from absorption column.)
    origin.Execute(f'plotxy (col(1), col({ycol})) plot:=200;')

    safe_title = title.replace('"', "'")
    origin.Execute(f'page.longname$ = "{safe_title}";')
    origin.Execute('layer.x.title$ = "Wavelength (nm)";')
    origin.Execute(f'layer.y.title$ = "{ylab}";')
    origin.Execute('legend.update;')


def save_project(origin, save_path):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    origin.Execute(f'save -n "{lt_escape_path(save_path)}";')


# =========================
# OVERLAY SUPPORT (resample to common x grid)
# =========================
def lin_interp(x, y, xq):
    """
    Simple linear interpolation for monotonic-ish x arrays.
    Returns float('nan') outside range.
    """
    if not x or not y or len(x) != len(y):
        return [float("nan")] * len(xq)

    # Ensure sorted by x (many are already sorted)
    pairs = [(xx, yy) for xx, yy in zip(x, y) if math.isfinite(xx) and math.isfinite(yy)]
    if len(pairs) < 5:
        return [float("nan")] * len(xq)

    pairs.sort(key=lambda p: p[0])
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    out = []
    j = 0
    for xv in xq:
        if xv < xs[0] or xv > xs[-1]:
            out.append(float("nan"))
            continue
        while j < len(xs) - 2 and xs[j+1] < xv:
            j += 1
        x0, x1 = xs[j], xs[j+1]
        y0, y1 = ys[j], ys[j+1]
        if abs(x1 - x0) < 1e-15:
            out.append(y0)
        else:
            t = (xv - x0) / (x1 - x0)
            out.append(y0 + t * (y1 - y0))
    return out


def build_overlay_grid(selected_blocks):
    """
    Grid based on overlapping region of all selected spectra in this batch.
    """
    mins = []
    maxs = []
    for b in selected_blocks:
        # choose any available spectrum for x-range
        spec = b.get("Absorption") or b.get("Emission")
        if not spec:
            continue
        x = spec["x"]
        if not x:
            continue
        mins.append(min(x))
        maxs.append(max(x))

    if not mins or not maxs:
        return None

    lo = max(mins)   # overlap start
    hi = min(maxs)   # overlap end
    if hi <= lo:
        # fallback: just use global min/max if no overlap
        lo = min(mins)
        hi = max(maxs)

    step = float(OVERLAY_GRID_STEP_NM)
    n = int((hi - lo) / step) + 1
    return [lo + i * step for i in range(n)]


def export_overlay_workbook(origin, overlay_name, proteins_in_batch, proteins_dict, which="Absorption"):
    """
    Make one workbook where:
      col1 = common wavelength grid
      col2.. = each protein's normalised spectrum (abs or em)
    """
    blocks = [proteins_dict[p] for p in proteins_in_batch]
    grid = build_overlay_grid(blocks)
    if not grid or len(grid) < 10:
        return None

    # Build table rows
    cols = [grid]
    labels = ["Wavelength_nm"]

    for p in proteins_in_batch:
        spec = proteins_dict[p].get(which)
        if not spec:
            continue
        y = spec["y_norm"]
        yq = lin_interp(spec["x"], y, grid)
        cols.append(yq)
        labels.append(p)

    if len(cols) < 2:
        return None

    # write TSV and import into Origin
    short = safe_page_name(overlay_name)
    temp_dir = os.path.join(find_project_root(), "Project csv data", "_origin_temp_export")
    tsv_path = os.path.join(temp_dir, f"{short}.tsv")

    with open(tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(labels)
        for i in range(len(grid)):
            row = []
            for c in cols:
                v = c[i]
                if isinstance(v, float) and (math.isnan(v) or not math.isfinite(v)):
                    row.append("")
                else:
                    row.append(str(v))
            w.writerow(row)

    # create workbook + import
    origin.Execute('newbook name:=TmpBook sheet:=1 option:=lsname;')
    origin.Execute(f'page.name$="{short}"; page.longname$="{overlay_name.replace(chr(34), chr(39))}";')
    origin.Execute(f'impASC fname:="{lt_escape_path(tsv_path)}" options:="TNAME:0,DELIM:9,ONAME:1";')

    # Plot all columns (2..N) vs col1
    ncols = len(labels)
    if ncols >= 3:
        origin.Execute(f'plotxy (col(1), col(2):col({ncols})) plot:=200;')
    else:
        origin.Execute('plotxy (col(1), col(2)) plot:=200;')

    origin.Execute('layer.x.title$="Wavelength (nm)";')
    origin.Execute('layer.y.title$="Normalised Intensity (0–1)";')
    origin.Execute('legend.update;')

    return short


# =========================
# MAIN
# =========================
def main():
    project_root = find_project_root()
    if not project_root:
        print("❌ Could not find project root automatically.")
        return

    master_csv = os.path.join(project_root, MASTER_REL_PATH)
    if not os.path.exists(master_csv):
        print("❌ Master CSV not found:", master_csv)
        return

    progress_path = os.path.join(project_root, PROGRESS_REL_PATH)

    proteins = load_spectra_from_master(master_csv)
    all_names = sorted(proteins.keys())

    prog = load_progress(progress_path)
    exported = set(prog.get("exported_proteins", []))

    remaining = [p for p in all_names if p not in exported]

    print("=== Master CSV → Origin exporter (no repeat, 2 graphs + overlays) ===")
    print("PROJECT_ROOT:", project_root)
    print("MASTER_CSV:", master_csv)
    print("Loaded proteins:", len(all_names))
    print("Already exported:", len(exported))
    print("Remaining:", len(remaining))

    if not remaining:
        print("✅ Nothing left to export. Delete progress file to restart:")
        print("  ", progress_path)
        return

    # pick next batch
    batch = remaining[:BATCH_SIZE] if BATCH_SIZE and BATCH_SIZE > 0 else remaining

    # output file name increments by how many batches we've already done
    batch_index = (len(exported) // max(1, BATCH_SIZE)) + 1
    out_name = f"{EXPORT_NAME_BASE}_batch{batch_index:02d}.opju"
    export_opju = os.path.join(project_root, out_name)

    temp_dir = os.path.join(project_root, "Project csv data", "_origin_temp_export")

    print("\nExporting proteins this run:", len(batch))
    print("Saving OPJU to:", export_opju)

    print("\nLaunching Origin...")
    origin = get_origin_instance(visible=True)
    origin.Execute("doc -s; doc -n;")

    total_books = 0
    total_specs = 0

    for protein in batch:
        block = proteins[protein]
        abs_spec = block.get("Absorption")
        em_spec = block.get("Emission")

        ref = abs_spec if abs_spec else em_spec
        if not ref:
            continue

        x = ref["x"]
        n = len(x)
        nan_col = [float("nan")] * n

        abs_raw = (abs_spec["y_raw"][:n] if abs_spec else nan_col)
        abs_nrm = (abs_spec["y_norm"][:n] if abs_spec else nan_col)
        em_raw  = (em_spec["y_raw"][:n] if em_spec else nan_col)
        em_nrm  = (em_spec["y_norm"][:n] if em_spec else nan_col)

        safe_short = safe_page_name(protein)
        tsv_path = os.path.join(temp_dir, f"{safe_short}.tsv")
        write_temp_tsv(tsv_path, x, abs_raw, abs_nrm, em_raw, em_nrm)

        book = create_workbook_and_import_tsv(origin, protein, tsv_path)

        # Create two correct graphs (only if the spectrum exists)
        if abs_spec is not None:
            create_single_graph(origin, book, protein, which="Absorption", plot_norm=PLOT_NORMALISED)
        if em_spec is not None:
            create_single_graph(origin, book, protein, which="Emission", plot_norm=PLOT_NORMALISED)

        total_books += 1
        total_specs += (1 if abs_spec else 0) + (1 if em_spec else 0)
        time.sleep(0.05)

    # Overlay workbooks/graphs for the batch (normalised)
    export_overlay_workbook(origin, f"Overlay_Abs_batch{batch_index:02d}", batch, proteins, which="Absorption")
    export_overlay_workbook(origin, f"Overlay_Em_batch{batch_index:02d}", batch, proteins, which="Emission")

    save_project(origin, export_opju)

    # update progress
    prog.setdefault("exported_proteins", [])
    prog["exported_proteins"].extend(batch)
    # unique-preserve order
    seen = set()
    new_list = []
    for p in prog["exported_proteins"]:
        if p not in seen:
            seen.add(p)
            new_list.append(p)
    prog["exported_proteins"] = new_list
    save_progress(progress_path, prog)

    print("\n=== Export complete ===")
    print("Exported proteins:", total_books)
    print("Exported spectra: ", total_specs)
    print("Saved OPJU to:\n ", export_opju)
    print("\nProgress file:\n ", progress_path)
    print("Run again to export the next batch (it will NOT repeat).")


if __name__ == "__main__":
    main()
