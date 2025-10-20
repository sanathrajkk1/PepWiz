from pathlib import Path
import re

def open_reader(ms_path: Path):
    """Open an mzML or mzXML file using pyteomics."""
    from pyteomics import mzxml, mzml
    opener = mzml.MzML if ms_path.suffix.lower() == ".mzml" else mzxml.MzXML
    return opener(str(ms_path))

def precursor_mz_from_spec(spec):
    """Extract precursor m/z from an MS2 spectrum (mzML or mzXML)."""
    try:
        if 'precursorList' in spec:  # mzML
            sel = spec['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
            return float(sel.get('selected ion m/z') or sel.get('selected ion m/z array', [None])[0])
        elif 'precursorMz' in spec:   # mzXML
            return float(spec['precursorMz'][0]['precursorMz'])
    except Exception:
        return None

def list_precursors_with_counts(ms_path: Path, dedup_ppm: float = 10.0):
    """List unique precursor m/z clusters and counts."""
    parents = []
    try:
        with open_reader(ms_path) as reader:
            for spec in reader:
                ms_level = spec.get('ms level') or spec.get('msLevel')
                try:
                    if int(ms_level) != 2:
                        continue
                except Exception:
                    continue
                pmz = precursor_mz_from_spec(spec)
                if pmz is not None:
                    parents.append(pmz)
    except ImportError as e:
        raise RuntimeError("pyteomics (and lxml) are required. Run:\n  py -m pip install pyteomics lxml") from e

    if not parents:
        return []
    parents.sort()
    clusters, cur = [], [parents[0]]
    for mz in parents[1:]:
        center = sum(cur)/len(cur)
        if abs(mz-center)/center*1e6 <= dedup_ppm:
            cur.append(mz)
        else:
            clusters.append({"mz": sum(cur)/len(cur), "count": len(cur)})
            cur = [mz]
    clusters.append({"mz": sum(cur)/len(cur), "count": len(cur)})
    clusters.sort(key=lambda r: (-r["count"], r["mz"]))
    return clusters

def iter_filtered_ms2_peaks(ms_path: Path, precursor_mz: float | None, ppm_tol: float,
                            rt_min: float | None, rt_max: float | None):
   
    with open_reader(ms_path) as reader:
        for spec in reader:
            ms_level = spec.get('ms level') or spec.get('msLevel')
            try:
                if int(ms_level) != 2:
                    continue
            except Exception:
                continue

            # RT filter (if available)
            if rt_min is not None or rt_max is not None:
                # mzML: 'scanList'/'scan'/'scan start time' (minutes). mzXML sometimes 'retentionTime' in seconds (PTxxS).
                rt = spec.get('scanList', {}).get('scan', [{}])[0].get('scan start time')
                if rt is None:
                    # try mzXML ISO8601 duration 'PTxxS'
                    iso = spec.get('retentionTime')
                    if isinstance(iso, str) and iso.startswith("PT") and iso.endswith("S"):
                        try:
                            rt = float(iso[2:-1]) / 60.0
                        except Exception:
                            rt = None
                # apply window if we have RT
                if rt is not None:
                    if rt_min is not None and rt < rt_min: 
                        continue
                    if rt_max is not None and rt > rt_max: 
                        continue

            # Precursor m/z filter (if provided)
            if precursor_mz is not None:
                pmz = precursor_mz_from_spec(spec)
                if pmz is None:
                    continue
                if abs(pmz - precursor_mz) / precursor_mz * 1e6 > ppm_tol:
                    continue

            mzs = spec['m/z array']; ints = spec['intensity array']
            yield list(zip(mzs, ints))

def average_spectrum(peaks_lists, bin_ppm: float = 10.0, top_n: int | None = 200):
    """Average multiple MS2 peak lists into one centroided spectrum."""
    if not peaks_lists:
        return []
    all_peaks = []
    for peaks in peaks_lists:
        all_peaks.extend(peaks)
    if not all_peaks:
        return []
    all_peaks.sort(key=lambda x: x[0])

    out = []
    i, N = 0, len(all_peaks)
    while i < N:
        mz0, I0 = all_peaks[i]
        lo, hi = mz0*(1-bin_ppm*1e-6), mz0*(1+bin_ppm*1e-6)
        j, sum_I, sum_mzI = i, 0.0, 0.0
        while j < N and lo <= all_peaks[j][0] <= hi:
            mzj, Ij = all_peaks[j]
            sum_I += Ij
            sum_mzI += mzj * Ij
            j += 1
        if sum_I > 0:
            out.append((sum_mzI/sum_I, sum_I))
        i = j

    if not out:
        return out
    if isinstance(top_n, int) and top_n > 0 and len(out) > top_n:
        out = sorted(out, key=lambda x: x[1], reverse=True)[:top_n]
    out.sort(key=lambda x: x[0])
    return out

