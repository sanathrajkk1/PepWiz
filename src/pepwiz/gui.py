#!/usr/bin/env python3
from __future__ import annotations

#1) Imports:

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
import re
import os, json, shutil, subprocess, tempfile

from pepwiz.msconvert_utils import (
    find_msconvert,
    run_msconvert,
)

from pepwiz.mzml_utils import (
    open_reader,
    precursor_mz_from_spec,
    list_precursors_with_counts,
    iter_filtered_ms2_peaks,
    average_spectrum,
)

from pepwiz.match_engine import (
    PROTON,
    WATER,
    AA_MASS,
    generate_theoretical_by,
    ppm_error,
    ppm_delta,
    calc_fragments,
    legacy_summary_from_spectrum,
    nearest_match,
    ion_meta,
    compute_cleavages_from_masses,
)

from pepwiz.visualize import (
    export_fragment_image,
    export_annotated_spectrum,
)

from pepwiz.io_legacy import write_legacy_out
  
       
def parse_rt_window(s: str):
    """Parse 'min,max' minutes string to (rt_min, rt_max) floats or (None, None)."""
    if not s or "," not in s:
        return None, None
    a, b = s.split(",", 1)
    try:
        rt_min = float(a.strip()) if a.strip() else None
    except ValueError:
        rt_min = None
    try:
        rt_max = float(b.strip()) if b.strip() else None
    except ValueError:
        rt_max = None
    return rt_min, rt_max


class PepWizGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PepWiz")         # window title text
        self.geometry("800x500")     #default size
      

        self.msfile_var = tk.StringVar()           #path to mzXML/mzML
        self.seq_var    = tk.StringVar()           #peptide suquence
        self.ppm_var    = tk.StringVar()           #ppm error
        self.z_var      = tk.StringVar()           #fragment charge
        self.precursor_var = tk.StringVar()        #MS1 precursor m/z you filtered on
        self.rt_var        = tk.StringVar()        #retention time
        self.B_var      = tk.StringVar()
        self.J_var      = tk.StringVar()
        self.X_var      = tk.StringVar()
        self.keep_mzml_var = tk.BooleanVar(value=False)  # default: don’t keep mzML
        self.term_mod_var = tk.StringVar(value="None")
        self.topn_var = tk.StringVar(value="200")
        self.draw_img_var = tk.BooleanVar(value=False)   #default dont generate fragment image
        
        self.ppm_var.set("10")
        self.draw_spec_var  = tk.BooleanVar(value=False)  # export annotated spectrum
        self.label_topn_var = tk.StringVar(value="30")    # max labels to draw on spectrum
        self.min_pct_var    = tk.StringVar(value="5")     # min % of base peak to label

#4) Widgets and layout
        pad = {"padx": 8, "pady": 6}
        
        #-------Top row: file chooser-------
        top = ttk.Frame(self); top.pack(fill=tk.X, **pad)
        ttk.Label(top, text="Choose .raw/.mzXML/.mzML file:").grid(row=0, column=0, sticky=tk.W)
        ttk.Entry(top, textvariable=self.msfile_var, width=70).grid(row=0, column=1, sticky=tk.W)
        ttk.Button(top, text="Browse", command=self._choose_msfile).grid(row=0, column=2, sticky=tk.W)
        
        keep_row = ttk.Frame(self); keep_row.pack(fill=tk.X, padx=8, pady=0)
        ttk.Checkbutton(
            keep_row,
            text="Save converted mzML next to RAW",
            variable=self.keep_mzml_var
        ).grid(row=0, column=0, sticky=tk.W)

        #-------Middle row: parameteres------
        mid = ttk.Frame(self); mid.pack(fill=tk.X, **pad)
        ttk.Label(mid, text="Peptide sequence:").grid(row=0, column=0, sticky=tk.W)
        ttk.Entry(mid, textvariable=self.seq_var, width=30).grid(row=0, column=1, sticky=tk.W)
        
        mod_row = ttk.Frame(self); mod_row.pack(fill=tk.X, padx=8, pady=0)
        ttk.Label(mod_row, text="Terminal modification:").grid(row=0, column=0, sticky=tk.W)
        ttk.OptionMenu(
            mod_row,
            self.term_mod_var,
            self.term_mod_var.get(),
            "None",
            "C-term: Amidated",
            "C-term: Dehydrated",
            "C-term: Decarboxylated (Daptides)"
        ).grid(row=0, column=1, sticky=tk.W)
        
        ttk.Label(mid, text="PPM tolerance:").grid(row=1, column=0, sticky=tk.W)
        ttk.Entry(mid, textvariable=self.ppm_var, width=10).grid(row=1, column=1, sticky=tk.W)

        ttk.Label(mid, text="Fragment charge(s):").grid(row=2, column=0, sticky=tk.W)
        ttk.Entry(mid, textvariable=self.z_var, width=10).grid(row=2, column=1, sticky=tk.W)
        ttk.Label(mid, text="e.g., 1 or 1,2").grid(row=2, column=2, sticky=tk.W)

        topn_row = ttk.Frame(self); topn_row.pack(fill=tk.X, padx=8, pady=0)
        ttk.Label(topn_row, text="Top peaks to match:").grid(row=0, column=0, sticky=tk.W)
        ttk.Entry(topn_row, textvariable=self.topn_var, width=8).grid(row=0, column=1, sticky=tk.W)

        img_row = ttk.Frame(self); img_row.pack(fill=tk.X, padx=8, pady=0)
        ttk.Checkbutton(img_row, text="Export fragment coverage image (optional)",
                        variable=self.draw_img_var).grid(row=0, column=0, sticky=tk.W)
        
        row_spec = ttk.Frame(self); row_spec.pack(fill=tk.X, padx=8, pady=0)
        ttk.Checkbutton(
            row_spec,
            text="Export annotated MS/MS spectrum",
            variable=self.draw_spec_var
        ).grid(row=0, column=0, sticky=tk.W)
        
        flt = ttk.Frame(self); flt.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(flt, text="Precursor m/z (most abundant MS1):").grid(row=0, column=0, sticky=tk.W)
        ttk.Entry(flt, textvariable=self.precursor_var, width=14).grid(row=0, column=1, sticky=tk.W)
        ttk.Label(flt, text="RT window min–max (min, optional):").grid(row=0, column=2, sticky=tk.W, padx=(16,0))
        ttk.Entry(flt, textvariable=self.rt_var, width=14).grid(row=0, column=3, sticky=tk.W)
        
        cust = ttk.Frame(self); cust.pack(fill=tk.X, padx=8, pady=6)
        ttk.Label(cust, text="Custom residue masses (Da, optional):").grid(row=0, column=0, columnspan=6, sticky=tk.W)
    
        ttk.Label(cust, text="B:").grid(row=1, column=0, sticky=tk.E)
        ttk.Entry(cust, textvariable=self.B_var, width=10).grid(row=1, column=1, sticky=tk.W)

        ttk.Label(cust, text="J:").grid(row=1, column=2, sticky=tk.E)
        ttk.Entry(cust, textvariable=self.J_var, width=10).grid(row=1, column=3, sticky=tk.W)

        ttk.Label(cust, text="X:").grid(row=1, column=4, sticky=tk.E)
        ttk.Entry(cust, textvariable=self.X_var, width=10).grid(row=1, column=5, sticky=tk.W)

        ttk.Separator(self).pack(fill=tk.X, **pad)
        
        #----Controls: Run + Progress
        ctrls = ttk.Frame(self); ctrls.pack(fill=tk.X, **pad)
        self.run_btn = ttk.Button(ctrls, text="Run", command=self._on_run)
        self.run_btn.pack(side=tk.LEFT)
        self.prog = ttk.Progressbar(ctrls, mode="determinate")
        self.prog.pack(fill=tk.X, expand=True, side=tk.LEFT, padx=8)
        
        ttk.Separator(self).pack(fill=tk.X, **pad)
        
       
         # ------- Log box -------
        logf = ttk.Frame(self); logf.pack(fill=tk.BOTH, expand=True, **pad)
        ttk.Label(logf, text="Log").pack(anchor="w")
        self.log_text = tk.Text(logf, height=12)
        self.log_text.pack(fill=tk.BOTH, expand=True)
        
        
        exe = find_msconvert(None)
        self._log(f"msconvert: {'found at ' + exe if exe else 'not found (will prompt on RAW)'}")
        self._log("Dependencies: install with `py -m pip install pyteomics lxml`")
        
    # Helper: file open dialog
    def _choose_msfile(self):
        p = filedialog.askopenfilename(
            title="Choose RAW/mzML/mzXML",
            filetypes=[
                ("Thermo RAW", "*.raw"),
                ("mzML", "*.mzML"),
                ("mzXML", "*.mzXML"),
                ("All files", "*.*"),
            ],
        )
        if not p:
            return

        src = Path(p)
        self._source_for_output = src  # remember original selection for .out location

        if src.suffix.lower() == ".raw":
            try:
                keep = bool(self.keep_mzml_var.get())
                dest = src.parent if keep else None
                mzml_path = run_msconvert(src, dest_dir=dest, log_fn=self._log, overwrite=True)
                self.msfile_var.set(str(mzml_path))
                target_for_listing = mzml_path
            except RuntimeError as e:
                self._log(str(e))
                messagebox.showerror("msconvert error", str(e))
                return
        else:
            self.msfile_var.set(str(src))
            target_for_listing = src


        # Print parent-ion list to the log (top 20)
        try:
            clusters = list_precursors_with_counts(target_for_listing, dedup_ppm=10.0)
            if not clusters:
                self._log("No MS2 parent ions found.")
            else:
                self._log(f"Found {sum(c['count'] for c in clusters)} MS2 scans grouped into {len(clusters)} parent ions:")
                header = f"{'Rank':>4}  {'Parent m/z':>12}  {'MS2 scans':>9}"
                self._log(header); self._log("-" * len(header))
                for i, c in enumerate(clusters[:20], 1):
                    self._log(f"{i:>4}  {c['mz']:>12.4f}  {c['count']:>9}")
                if len(clusters) > 20:
                    self._log(f"... and {len(clusters) - 20} more")
                self._log("Tip: copy a parent m/z into the 'Precursor m/z' box above.")
        except RuntimeError as e:
            self._log(str(e))


    # Helper: append to log
    def _log(self, msg: str):
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)

        # Run button handler
    def _on_run(self):
        self.prog.config(mode="indeterminate"); self.prog.start(12)
        try:
            msfile = Path(self.msfile_var.get()).expanduser()
            if not msfile.exists():
                messagebox.showerror("Missing file", "Choose a valid RAW/mzML/mzXML file.")
                return

            # If RAW, convert before analysis. Save next to RAW if checkbox is on.
            if msfile.suffix.lower() == ".raw":
                try:
                    keep = bool(self.keep_mzml_var.get())
                    dest = msfile.parent if keep else None
                    msfile = run_msconvert(msfile, dest_dir=dest, log_fn=self._log)
                except RuntimeError:
                    self._log("msconvert not found or failed to run.\n")
                    self._log("PepWiz requires ProteoWizard's msconvert.exe to process RAW files.\n")

                    self._log("To locate msconvert.exe, open PowerShell and run:")
                    self._log("  Get-ChildItem -Path 'C:\\' -Filter msconvert.exe -Recurse -ErrorAction SilentlyContinue\n")

                    self._log("Once you find the full path (for example: C:\\Program Files\\ProteoWizard\\msconvert.exe), set it permanently by running:")
                    self._log("  setx PEPWIZ_MS_CONVERT \"C:\\\\Program Files\\\\ProteoWizard\\\\msconvert.exe\"\n")

                    self._log("After setting it, restart PepWiz — it will automatically use the path from the environment variable.")
                    self._log("If you prefer, you can also just use mzML or mzXML files directly instead of RAW.\n")

                    messagebox.showerror(
                        "msconvert missing",
                        "PepWiz could not locate ProteoWizard msconvert.exe.\n"
                        "Check the log below for step-by-step PowerShell commands to fix it."
                    )
                    return



            seq = self.seq_var.get().strip().upper()
            if not seq:
                messagebox.showerror("Missing sequence", "Please enter a peptide sequence.")
                return

            try:
                ppm = float(self.ppm_var.get())
            except Exception:
                messagebox.showerror("Invalid PPM", "PPM tolerance must be a number.")
                return
            if ppm <= 0:
                messagebox.showerror("Invalid PPM", "PPM tolerance must be > 0.")
                return
                
            try:
                charges = [int(z.strip()) for z in self.z_var.get().split(",") if z.strip()]
                if not charges:
                    raise ValueError
            except Exception:
                messagebox.showerror("Invalid charges", "Use integers like 1 or 1,2")
                return
            if any(z <= 0 for z in charges):
                messagebox.showerror("Invalid charges", "Fragment charges must be positive integers.")
                return

            # --- Optional B/J/X overrides ---
            overrides = {}
            def _maybe_float(s: str):
                s = s.strip()
                return float(s) if s else None

            try:
                b = _maybe_float(self.B_var.get())
                j = _maybe_float(self.J_var.get())
                x = _maybe_float(self.X_var.get())
            except ValueError:
                messagebox.showerror("Invalid mass", "B/J/X must be numbers.")
                return

            if b is not None: overrides['B'] = b
            if j is not None: overrides['J'] = j
            if x is not None: overrides['X'] = x

            # Guard: if sequence uses B/J/X but mass missing
            for letter in ('B','J','X'):
                if letter in seq and letter not in overrides:
                    messagebox.showerror("Missing mass", f"You used '{letter}' in the sequence but left its mass blank.")
                    return

            # Guard: if sequence uses anything other than B/J/X 
            unknown = sorted({aa for aa in seq if aa not in AA_MASS})
            unknown = [u for u in unknown if u not in overrides]  # B/J/X allowed if overridden
            if unknown:
                messagebox.showerror("Unknown residue(s)", f"Unknown residue(s): {', '.join(unknown)}")
                return

            self._log("Inputs look good!!")
            self._log(f"File: {msfile}")
            self._log(f"Sequence: {seq}")
            self._log(f"PPM: {ppm}, Charges: {charges}")

            # --- Build theoretical fragments ---
            mod_choice = self.term_mod_var.get()
            theo = calc_fragments(seq, charges, overrides, mod_choice)
            self._log(f"Terminal modification: {mod_choice}")



            # --- Optional precursor and RT filters from UI ---
            precursor_str = self.precursor_var.get().strip()
            rt_min, rt_max = parse_rt_window(self.rt_var.get())
            
            # Require and parse the user-provided precursor m/z (as the gate hint)
            if not precursor_str:
                messagebox.showerror("Missing precursor m/z", "Enter the precursor m/z you selected in Xcalibur.")
                return
            try:
                precursor_target = float(precursor_str)
            except ValueError:
                messagebox.showerror("Invalid precursor m/z", "Enter a number like 678.3456")
                return
            self._log(f"Precursor input: '{precursor_str}' -> {precursor_target:.4f}")
            # Build/refresh the parent-ion inventory and snap to nearest cluster center
            try:
                clusters = list_precursors_with_counts(msfile, dedup_ppm=10.0)
            except Exception:
                clusters = []

            snapped_parent = precursor_target
            snap_note = ""
            if clusters:
                nearest = min(clusters, key=lambda c: ppm_delta(precursor_target, c["mz"]))
                delta_ppm = ppm_delta(precursor_target, nearest["mz"])
                snapped_parent = nearest["mz"]
                snap_note = f"(snapped to {snapped_parent:.4f}, Δ={delta_ppm:.2f} ppm, scans={nearest['count']})"
                if delta_ppm > 50:
                    self._log(f"Warning: entered parent {precursor_target:.4f} is {delta_ppm:.1f} ppm from nearest cluster {snapped_parent:.4f}.")
            else:
                self._log("No parent clusters found; using the typed precursor m/z as-is.")

            self._log(f"Precursor gate centered at {snapped_parent:.4f} {snap_note}")

            # Collect filtered MS2 scans (use the *float* snapped_parent)
            filtered = list(iter_filtered_ms2_peaks(msfile, snapped_parent, ppm, rt_min, rt_max))
            if not filtered:
                self._log("No MS2 scans passed the filters.")
                messagebox.showwarning("No scans", "No MS2 scans matched the precursor/RT filters.")
                return

            try:
                top_n = int(self.topn_var.get().strip() or "200")
                if top_n <= 0:
                    top_n = 200
            except Exception:
                top_n = 200


            # --- Average them to one spectrum ---
            avg_spec = average_spectrum(filtered, bin_ppm=ppm, top_n=top_n)


            # --- Match against the averaged spectrum ---
            summary_rows = legacy_summary_from_spectrum(avg_spec, theo, ppm)

            # --- Write compact .out ---
           # Choose the base path for output: original selection if available, else the working file
            base_for_out = getattr(self, "_source_for_output", msfile)

            z_label = ",".join(str(z) for z in charges)
            out_base = base_for_out.with_suffix("")               # strip .raw/.mzML/.mzXML
            out_path = out_base.parent / f"{out_base.name}.z{z_label}.out"
            img_base = base_for_out.with_suffix("")  
            z_label  = ",".join(str(z) for z in charges)


            parent_mz     = snapped_parent
            scans_count   = len(filtered)
            bin_width_ppm = ppm

            # If you also want to embed the inventory in the .out:
            parent_clusters = None
            try:
                parent_clusters = list_precursors_with_counts(msfile, dedup_ppm=10.0)
            except Exception:
                pass

            write_legacy_out(
                out_path, seq, charges, summary_rows,
                parent_mz=snapped_parent,
                ppm_gate=ppm,
                scans_count=len(filtered),
                rt_min=rt_min, rt_max=rt_max,
                bin_ppm=ppm,
                parent_clusters=parent_clusters,
                term_mod=mod_choice,
            )


            self._log(f"Wrote legacy summary:\n{out_path}")
            
                     
            try:
                label_top_n = int(self.label_topn_var.get().strip() or "30")
            except Exception:
                label_top_n = 30
            try:
                min_pct = float(self.min_pct_var.get().strip() or "5")
            except Exception:
                min_pct = 5.0

            # Fragment map (SVG only)
            if self.draw_img_var.get() and summary_rows:
                img_base = base_for_out.with_suffix("")  # align with .out naming
                z_label = ",".join(str(z) for z in charges)
                frag_svg = img_base.parent / f"{img_base.name}.z{z_label}.fragments.svg"
                export_fragment_image(seq, summary_rows, frag_svg, log_fn=self._log)

            # Annotated spectrum (SVG only)
            if self.draw_spec_var.get() and avg_spec and summary_rows:
                img_base = base_for_out.with_suffix("")
                z_label = ",".join(str(z) for z in charges)
                spec_svg = img_base.parent / f"{img_base.name}.z{z_label}.spectrum.svg"
                export_annotated_spectrum(
                    avg_spec,
                    summary_rows,
                    spec_svg,
                    label_top_n=label_top_n,
                    min_pct=min_pct,
                    log_fn=self._log
                )

            
            else:
                if not self.draw_spec_var.get():
                    self._log("Spectrum export skipped: toggle is OFF.")
                elif not avg_spec:
                    self._log("Spectrum export skipped: no averaged spectrum (check filters/precursor).")
                elif not summary_rows:
                    self._log("Spectrum export skipped: no matched fragments to annotate.")

            
            # If we wrote mzML to a temp folder AND the user didn't ask to keep it, clean up
            try:
                if not self.keep_mzml_var.get():
                    tmp_root = Path(msfile).parent
                    if str(tmp_root).startswith(str(Path(tempfile.gettempdir()))):
                        shutil.rmtree(tmp_root, ignore_errors=True)
                        self._log("Cleaned up temporary mzML folder.")
            except Exception:
                pass
        finally:
             self.prog.stop()
             self.prog.config(mode="determinate", value=0)


# Standard bootstrap
if __name__ == "__main__":
    app = PepWizGUI()
    app.mainloop()
       