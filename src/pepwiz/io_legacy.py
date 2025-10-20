from __future__ import annotations
from pathlib import Path

def write_legacy_out(
    out_path: Path,
    peptide: str,
    charges,
    rows,                         # [{z, itype, idx, ion, theo, obs, ppm}]
    parent_mz: float | None = None,
    ppm_gate: float | None = None,
    scans_count: int | None = None,
    rt_min: float | None = None,
    rt_max: float | None = None,
    bin_ppm: float | None = None,
    parent_clusters=None,         # optional list from list_precursors_with_counts()
    term_mod: str | None = None,  # <--- ADD THIS
):
    z_label = ",".join(str(z) for z in charges)
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        # Header
        fh.write(f"{peptide} at fragment charge(s) z = {z_label}\n")

        # Run parameters / provenance
        details = []
        if parent_mz is not None: details.append(f"Parent m/z = {parent_mz:.4f}")
        if ppm_gate is not None:  details.append(f"precursor gate ±{ppm_gate} ppm")
        if scans_count is not None: details.append(f"scans averaged = {scans_count}")
        if rt_min is not None or rt_max is not None:
            details.append(f"RT window (min) = {'' if rt_min is None else f'{rt_min:.2f}'}-{'' if rt_max is None else f'{rt_max:.2f}'}")
        if bin_ppm is not None:   details.append(f"averaging bin = ±{bin_ppm} ppm")
        if term_mod and term_mod != "None":  # <--- include mod in header
            details.append(f"terminal mod = {term_mod}")
        details.append("grouping = by fragment charge (z), then b->y, then index")

        fh.write("Run parameters: " + " | ".join(details) + "\n")

        # Optional: parent-ion inventory...
        if parent_clusters:
            fh.write("Parent ions present (top 10 by MS2 count):\n")
            fh.write(f"{'Rank':>4}  {'Parent m/z':>12}  {'MS2 scans':>9}\n")
            fh.write("-" * 32 + "\n")
            for i, c in enumerate(parent_clusters[:10], 1):
                fh.write(f"{i:>4}  {c['mz']:>12.4f}  {c['count']:>9}\n")

        # Sections grouped by fragment charge
        frag_zs = sorted({r["z"] for r in rows if r["z"]})
        for fz in frag_zs:
            fh.write("\n")
            fh.write(f"[ Fragment charge z = {fz} ]\n")
            fh.write(f"{'IonType':6} {'Ion':10} {'Parent_m/z':>12} {'Pred_m/z':>12} {'Obs_m/z':>12} {'PPM':>8}\n")
            fh.write(f"{'-'*6} {'-'*10} {'-'*12} {'-'*12} {'-'*12} {'-'*8}\n")
            for r in rows:
                if r["z"] != fz:
                    continue
                itype = r["itype"]
                ion   = r["ion"]
                theo  = f"{r['theo']:.4f}"
                parent_col = f"{parent_mz:.4f}" if parent_mz is not None else "NA"
                if r["obs"] is None or r["ppm"] is None:
                    obs = "NA"; ppmv = "NA"
                else:
                    obs = f"{r['obs']:.4f}"; ppmv = f"{r['ppm']:.2f}"
                fh.write(f"{itype:6} {ion:10} {parent_col:>12} {theo:>12} {obs:>12} {ppmv:>8}\n")
