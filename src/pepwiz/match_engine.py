from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Tuple, Dict
import re

PROTON = 1.007276466812  # 
WATER  = 18.010564684    # 


H_ATOM             = 1.00784           
DECARB_DAPTIDE_NEU = 29.0022           # y-ion decarboxylation net (if that toggle is used)  
AMIDATION_DELTA    = -0.984016         # C-term amidation
DEHYDRATION_DELTA  = -18.010564684     # net -H2O 

# ---- Monoisotopic AA masses  ----
AA_MASS = {
    "A": 71.03711,  "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276,  "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
}

def _prefix_masses(seq: str) -> List[float]:
    acc = 0.0
    out: List[float] = []
    for aa in seq:
        acc += AA_MASS.get(aa.upper(), 0.0)
        out.append(acc)
    return out

def _suffix_masses(seq: str) -> List[float]:
    acc = 0.0
    out: List[float] = []
    for aa in reversed(seq):
        acc += AA_MASS.get(aa.upper(), 0.0)
        out.append(acc)
    return list(reversed(out))

def _mz(neutral_mass: float, z: int) -> float:
    return (neutral_mass + z * PROTON) / z

def _apply_terminal_mod(neutral_mass: float, term_mod: str | None, series: str) -> float:
    """
    Apply terminal modification to the neutral mass of a fragment if requested.
    term_mod values are the same strings your GUI uses (“None”, “Amidation”, “Dehydration”, “Daptide decarb (y)”, …).
    Only applies when it makes sense (e.g., y-series for decarb).
    """
    if not term_mod or term_mod == "None":
        return neutral_mass
    t = term_mod.lower().strip()
    if t.startswith("amid"):       # amidation (C-term)
        
        return neutral_mass + AMIDATION_DELTA
    if "dehydrat" in t:            # -H2O
        return neutral_mass + DEHYDRATION_DELTA
    if "decarb" in t and series == "y":  
        return neutral_mass - DECARB_DAPTIDE_NEU
    return neutral_mass

def generate_theoretical_by(
    seq: str,
    charges: Iterable[int],
    term_mod: str | None = None,
) -> List[Tuple[str, float]]:
    """
    Generate b/y series for the given fragment charge set.
    Returns [(ion_label, theo_mz)] like ('b5^2+', 345.1234).
    Mirrors your v12 “theoretical ion builder” behavior for b/y with an optional terminal mod applied. 
    """
    L = len(seq)
    if L < 2:
        return []
    pref = _prefix_masses(seq)          # b_n base
    suff = _suffix_masses(seq)          # y_n base (needs +H2O)

    theo: List[Tuple[str, float]] = []
    chs = sorted(int(z) for z in charges if int(z) > 0)
    if not chs:
        chs = [1]

    for i in range(1, L):               # cut between i and i+1
        # Neutral fragment masses (peptide fragments without protons)
        b_mass = pref[i-1]              # b_i
        y_mass = suff[i] + WATER        # y_{L-i} + H2O

        # Optional terminal adjustments
        b_mass_adj = _apply_terminal_mod(b_mass, term_mod, "b")
        y_mass_adj = _apply_terminal_mod(y_mass, term_mod, "y")

        for z in chs:
            theo.append((f"b{i}^{z}+", _mz(b_mass_adj, z)))
            y_idx = L - i
            theo.append((f"y{y_idx}^{z}+", _mz(y_mass_adj, z)))
    return theo
    
def ppm_error(obs_mz: float, theo_mz: float, *, signed: bool = True) -> float:
    """
    Return the PPM error between observed and theoretical.
    signed=True -> (obs - theo)/theo * 1e6; signed=False -> abs value.
    """
    if theo_mz == 0:
        return 0.0
    ppm = (obs_mz - theo_mz) / theo_mz * 1e6
    return ppm if signed else abs(ppm)

def ppm_delta(a: float, b: float) -> float:
    """Absolute ppm difference between two m/z values."""
    if b == 0:
        return 0.0
    return abs((a - b) / b) * 1e6

def calc_fragments(seq: str, charges, overrides: dict, term_mod_choice: str):
    # Residue masses, honoring overrides (e.g., B/J/X if the user set them)
    
    masses = [(overrides[aa] if aa in overrides else AA_MASS[aa]) for aa in seq]

    # b-prefix neutral masses
    b_prefix = []
    s = 0.0
    for i, m in enumerate(masses, start=1):
        s += m
        b_prefix.append((i, s))

    # y-suffix neutral masses
    y_suffix = []
    s = 0.0
    for i, m in enumerate(reversed(masses), start=1):
        s += m
        y_suffix.append((i, s))

    out = []
    for z in charges:
        z = int(z)

        # b-ions are unaffected by C-term mods
        for i, s in b_prefix:
            out.append((f"b{i}^{z}+", (s + z*PROTON)/z))

        # y-ions depend on terminal modification
        if term_mod_choice == "C-term: Amidated":
            for i, s in y_suffix:
                out.append((f"y{i}^{z}+", (s + WATER + AMIDATION_DELTA + z*PROTON)/z))

        elif term_mod_choice == "C-term: Dehydrated":
            for i, s in y_suffix:
                out.append((f"y{i}^{z}+", (s + WATER + DEHYDRATION_DELTA + z*PROTON)/z))

        elif term_mod_choice == "C-term: Decarboxylated (Daptides)":
            for i, s in y_suffix:
                out.append((f"y{i}^{z}+", (s - DECARB_DAPTIDE_NEU + H_ATOM + z*PROTON)/z))

        else:  # "None"
            for i, s in y_suffix:
                out.append((f"y{i}^{z}+", (s + WATER + z*PROTON)/z))

    return out
    
def ion_meta(ion_label: str):
    """
    Return (itype, idx, z) from labels like 'b5^2+' or 'y10^1+'.
    itype: 'b' or 'y'; idx: int; z: int (fragment charge)
    """
    m = re.match(r'^([by])(\d+)\^(\d+)\+$', ion_label)
    if not m:
        return ("?", 0, 0)
    return (m.group(1), int(m.group(2)), int(m.group(3)))

def nearest_match(spectrum: List[Tuple[float, float]], target_mz: float, ppm_tol: float):
    """
    Return (obs_mz, intensity, ppm_error) for the best hit within ±ppm_tol, else None.
    Mirrors your v12 matcher.  
    """
    if not spectrum:
        return None
    best = None
    best_ppm = float("inf")
    for mz, inten in spectrum:
        ppm = abs(mz - target_mz) / target_mz * 1e6
        if ppm <= ppm_tol and ppm < best_ppm:
            best_ppm = ppm
            best = (mz, inten, ppm)
    return best

def legacy_summary_from_spectrum(
    spectrum: List[Tuple[float, float]],
    theo_ions: List[Tuple[str, float]],
    ppm_tol: float
) -> List[Dict]:
    """
    Same rows/ordering as your v12 summary.  
    Returns: [{z, itype, idx, ion, theo, obs, ppm, inten}] sorted by z -> (b then y) -> idx.
    """
    rows: List[Dict] = []
    for ion_label, theo_mz in theo_ions:
        hit = nearest_match(spectrum, theo_mz, ppm_tol)
        itype, idx, z = ion_meta(ion_label)
        if hit is None:
            continue
        obs_mz, inten, ppm = hit
        rows.append({
            "z": z, "itype": itype, "idx": idx, "ion": ion_label,
            "theo": theo_mz, "obs": obs_mz, "ppm": ppm, "inten": inten
        })
    rows.sort(key=lambda r: (r["z"], 0 if r["itype"] == "b" else 1, r["idx"]))
    return rows

def compute_cleavages_from_masses(seq: str, matched_rows):
    """
    Return two sets of cleavage indices (between 1..len(seq)-1):
      b_cuts: positions i where b_i observed
      y_cuts: positions i where y_i observed (y_i == cut between len-i and len-i+1)
    """
    L = len(seq)
    b_cuts = set()
    y_cuts = set()
    for r in matched_rows:
        itype, idx, z = r["itype"], r["idx"], r["z"]
        if itype == "b" and 1 <= idx < L:
            b_cuts.add(idx)
        elif itype == "y" and 1 <= idx < L:
            # y_i corresponds to cut between L-i and L-i+1
            cut = L - idx
            if 1 <= cut < L:
                y_cuts.add(cut)
    return b_cuts, y_cuts