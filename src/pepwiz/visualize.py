from __future__ import annotations

def _ensure_matplotlib(log_fn=print):
 
    try:
        import importlib, sys
        m = importlib.import_module("matplotlib")
       
        want = {
            "svg.fonttype": "none",  
            "svg.hashsalt": "",       
            "svg.image_inline": True, 
            
        }
        for k, v in want.items():
            if k in m.rcParams:
                m.rcParams[k] = v

        m.use("Agg")  # headless backend for GUI-less savefig
        plt = importlib.import_module("matplotlib.pyplot")
        log_fn(f"matplotlib OK: v{getattr(m, '__version__', '?')} | exe={sys.executable}")
        return True, plt
    except Exception as e:
        import traceback, sys
        tb = traceback.format_exc(limit=3)
        log_fn(f"matplotlib import failed: {type(e).__name__}: {e}\n{tb}\nexe={sys.executable}")
        return False, None
     
        
def export_fragment_image(seq, matched_rows, out_svg, log_fn=print):
  
    ok, plt = _ensure_matplotlib(log_fn)
    if not ok:
        return

    import matplotlib as mpl
    from math import floor

    L = len(seq)
    if L < 2:
        log_fn("Fragment image: sequence too short to draw cut markers.")
        return

    # ==== layout params (tweak to taste) ====
    fs_seq  = max(14, min(22, 12 + floor(L / 6)))  # bold sequence font size
    fs_ion  = 12                                    # ion label font size
    v_half  = 0.85                                  # half-height of the mid-line
    h_stub  = 0.60                                  # horizontal stub length
    pad_txt = 0.05                                  # tiny offset for label from the stub
    linecol = "#6f8ea8"                             # neutral teal-grey for lines
    b_col, y_col = "red", "blue"                    # label colors

    # dashed-line pattern (dotted look)
    DOT = (0, (2, 3))  # (offset, on/off pattern) -> small dots

    # ==== figure ====
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * L), 2.6), dpi=150)

    # 1) peptide sequence (bold monospace)
    for i, aa in enumerate(seq, 1):
        ax.text(i, 0.0, aa,
                ha="center", va="center",
                fontsize=fs_seq, fontfamily="monospace", fontweight="bold")

    # helper: single mathtext label with superscript charge
    def ion_label(itype, idx, z):
        # show '+' for z=1, or 'z+' otherwise
        sup = "+" if (z is None or z == 1) else f"{z}+"
        return rf"${itype}_{{{idx}}}^{{{sup}}}$"

    # 2) draw matched ions
    
    for r in matched_rows:
        it, idx, z = r.get("itype"), r.get("idx"), r.get("z")
        if it not in ("b", "y") or not idx or not (1 <= idx < L):
            continue

        # cut position: b_n cut lies after residue n; y_n cut is after residue (L-n)
        x_cut = idx + 0.5 if it == "b" else (L - idx) + 0.5

        if it == "y":
            # dotted mid-line: only upwards from center
            ax.plot([x_cut, x_cut], [0.0, v_half], linestyle=DOT, color=linecol, lw=1.2)
            # dotted horizontal stub to the RIGHT
            x0, x1, y = x_cut, x_cut + h_stub, v_half
            ax.plot([x0, x1], [y, y], linestyle=DOT, color=linecol, lw=1.2)
            # label centered on the stub
            ax.text((x0 + x1) / 2.0, y + pad_txt, ion_label("y", idx, z),
                    color=y_col, ha="center", va="bottom", fontsize=fs_ion)
        else:
            # dotted mid-line: only downwards from center
            ax.plot([x_cut, x_cut], [0.0, -v_half], linestyle=DOT, color=linecol, lw=1.2)
            # dotted horizontal stub to the LEFT
            x0, x1, y = x_cut - h_stub, x_cut, -v_half
            ax.plot([x0, x1], [y, y], linestyle=DOT, color=linecol, lw=1.2)
            # label centered on the stub
            ax.text((x0 + x1) / 2.0, y - pad_txt, ion_label("b", idx, z),
                    color=b_col, ha="center", va="top", fontsize=fs_ion)

    # tidy canvas
    ax.axis("off")
    ax.set_xlim(0.5, L + 0.5)
    ax.set_ylim(-1.2, 1.2)

    fig.tight_layout()
    try:
        fig.savefig(out_svg, format="svg", bbox_inches="tight", transparent=True)
        log_fn(f"Saved fragment image SVG (editable): {out_svg}")
    except Exception as e:
        log_fn(f"Failed to save fragment SVG: {e}")
    finally:
        plt.close(fig)


def export_annotated_spectrum(
    avg_spec,
    matched_rows,
    out_svg,
    label_top_n=9999,   # kept for signature compat
    min_pct=0.0,        # unused for matched; compat
    log_fn=print
):
    ok, plt = _ensure_matplotlib(log_fn)
    if not ok:
        return

    try:
        if not avg_spec:
            log_fn("Spectrum export: no averaged spectrum to plot.")
            return

        from matplotlib.ticker import MultipleLocator

        # --- data prep ---
        mzs  = [p[0] for p in avg_spec]
        ints = [p[1] for p in avg_spec]
        base = max(ints) if ints else 1.0
        norm = [i / base * 100.0 for i in ints]

        def ppm(a, b):  # a=obs, b=ref
            return abs(a - b) / b * 1e6

        tol_ppm = 10.0  # for matching observed m/z to a bar in the averaged spectrum

        # Fast lookup for heights near an observed m/z
        def height_at(obs_mz, fallback_inten):
            for x, y in zip(mzs, norm):
                if ppm(x, obs_mz) <= tol_ppm:
                    return y
            return (fallback_inten / base * 100.0)

        # Collect matched ion m/z and per-ion color
        matched = []
        matched_mzs = []
        for r in matched_rows:
            obs = r.get("obs")
            it = r.get("itype")
            idx = r.get("idx")
            z   = r.get("z")
            if obs is None or it not in ("b", "y") or not idx or not z:
                continue
            matched.append((obs, it, idx, z, r.get("inten", 0.0)))
            matched_mzs.append(obs)

        # --- figure size ---
        span = (max(mzs) - min(mzs)) if mzs else 0.0
        width = max(8.0, min(16.0, span / 150.0 + 8.0))
        fig, ax = plt.subplots(figsize=(width, 4.2), dpi=150)

        # 1) draw ALL bars as neutral gray first
        for x, y in zip(mzs, norm):
            ax.vlines(x, 0.0, y, color="#A0A0A0", linewidth=0.9)

        # 2) overlay matched peaks with ion-specific colors
        #    (choose the bar height from the averaged spectrum)
        for obs, it, idx, z, inten in matched:
            h = height_at(obs, inten)
            col = "red" if it == "b" else "blue"
            ax.vlines(obs, 0.0, h, color=col, linewidth=1.2)

        # 3) annotate ALL matched ions with label (colored) + m/z (black)
        #    font scales with number of labels for readability
        n_labels = max(1, len(matched))
        ion_fs   = max(8, min(11, int(11 - 0.004 * n_labels)))   # clamp 8–11
        mz_fs    = max(7, ion_fs - 1)

        placed_label_x = []
        for obs, it, idx, z, inten in sorted(matched, key=lambda t: t[0]):
            h  = height_at(obs, inten)
            y0 = h + 2.0
            if any(abs(obs - prev) < 5 for prev in placed_label_x):
                y0 += 6.0
            placed_label_x.append(obs)

            lbl = f"${it}_{{{idx}}}^{{{z}+}}$"
            col = "red" if it == "b" else "blue"

            # colored ion label
            ax.text(obs, y0, lbl, color=col, ha="center", va="bottom", fontsize=ion_fs)
            # dotted black connector
            ax.vlines(obs, h, y0 - 0.5, color="black", linestyle=":", linewidth=0.9)
            # black m/z above
            ax.text(obs, y0 + 6.0, f"{obs:.4f}", color="black", ha="center", va="bottom", fontsize=mz_fs)

        # --- axes & limits ---
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")

        # X-range: focus on matched region (−50, +100); fallbacks sensible
        if matched_mzs:
            lo, hi = min(matched_mzs), max(matched_mzs)
            x_min, x_max = lo - 50.0, hi + 100.0
            if mzs:
                x_min = max(min(mzs), x_min)
                x_max = min(max(mzs), x_max)
        else:
            # fallback: show whole spectrum ±50/100 if we somehow have no matches
            if mzs:
                lo, hi = min(mzs), max(mzs)
                x_min, x_max = lo - 50.0, hi + 100.0
            else:
                x_min, x_max = 0.0, 1000.0

        # Ensure non-zero width
        if x_max - x_min < 50:
            mid = (x_min + x_max) / 2.0
            x_min, x_max = mid - 50.0, mid + 50.0

        ax.set_xlim(x_min, x_max)

        # Major ticks every 100 m/z
        ax.xaxis.set_major_locator(MultipleLocator(100))

        # Y: fixed 0–100%
        ax.set_ylim(0.0, 100.0)

        # Clean frame
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        fig.savefig(out_svg, format="svg", bbox_inches="tight", transparent=True)
        log_fn(f"Saved annotated spectrum SVG (editable): {out_svg}")
    except Exception as e:
        log_fn(f"Spectrum export failed: {type(e).__name__}: {e}")
    finally:
        try:
            plt.close(fig)
        except Exception:
            pass