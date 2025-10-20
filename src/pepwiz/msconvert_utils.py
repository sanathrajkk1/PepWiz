# msconvert_utils.py

from __future__ import annotations
import subprocess, os
from pathlib import Path

def _looks_like_gui_or_cache(p: str) -> bool:
    pl = p.lower()
    return (
        pl.endswith("msconvertgui.exe")
        or pl.endswith("msconvertgui_icon.exe")
        or "microsoft\\installer" in pl
    )

def _normalize_msconvert_path(p: str | None) -> str | None:
    if not p:
        return None
    pp = p
    # resolve .lnk if you kept _resolve_windows_shortcut; otherwise skip
    if pp.lower().endswith(".lnk"):
        try:
            pp = _resolve_windows_shortcut(pp) or ""
        except Exception:
            pp = ""
    if not pp:
        return None
    if _looks_like_gui_or_cache(pp):
        # try sibling/parent CLI
        sib = os.path.join(os.path.dirname(pp), "msconvert.exe")
        par = os.path.abspath(os.path.join(os.path.dirname(pp), "..", "msconvert.exe"))
        if os.path.exists(sib): return sib
        if os.path.exists(par): return par
        return None
    return pp if os.path.exists(pp) else None

def _preflight_msconvert(path: str) -> tuple[bool, str]:
    try:
        cp = subprocess.run([path, "--help"], capture_output=True, text=True, timeout=10)
        # ProteoWizard sometimes returns 0 or 1 for --help; both are fine if it printed usage.
        if cp.returncode in (0, 1):
            return True, ""
        return False, f"msconvert exited with code {cp.returncode}\nSTDOUT:\n{cp.stdout}\nSTDERR:\n{cp.stderr}"
    except Exception as e:
        return False, f"Failed to execute msconvert: {e}"

def find_msconvert(explicit_path: str | None = None) -> str | None:
    
    cand = _normalize_msconvert_path(explicit_path)
    if cand: return cand
    
    env = os.environ.get("PEPWIZ_MS_CONVERT") or os.environ.get("MSCONVERT_EXE")
    cand = _normalize_msconvert_path(env)
    if cand: return cand
  
    from shutil import which
    cand = _normalize_msconvert_path(which("msconvert.exe") or which("msconvert"))
    if cand: return cand
    
    default = r"C:\Program Files\ProteoWizard\msconvert.exe"
    cand = _normalize_msconvert_path(default)
    if cand: return cand

    return None

def run_msconvert(raw_path: str | Path, out_dir: str | Path | None = None, *,
                  overwrite: bool = False, extra_filters: list[str] | None = None,
                  log_fn=None, dest_dir: str | Path | None = None) -> Path:
    exe = find_msconvert(None)
    if not exe:
        msg = (
            "ProteoWizard `msconvert.exe` not found.\n\n"
            "To fix this:\n"
            " 1️. Install ProteoWizard:\n"
            "     https://proteowizard.sourceforge.io/download.html\n\n"
            " 2️. Locate the actual CLI binary (open Command Prompt and run):\n"
            "     where /r C:\\ msconvert.exe\n\n"
            "   • Ignore results like MSConvertGUI_Icon.exe — those are shortcuts.\n"
            " 3️. Register it permanently (Command Prompt):\n"
            "     setx PEPWIZ_MS_CONVERT \"<path to msconvert.exe>\"\n"
            " 4️. Restart PepWiz and try again.\n\n"
            "Tip: Verify the path prints usage (not a GUI):\n"
            "    \"%PEPWIZ_MS_CONVERT%\" --help\n"
        )
        if log_fn:
            log_fn(msg)
        raise RuntimeError(msg)

    ok, preflight_msg = _preflight_msconvert(exe)
    if not ok:
        err = (
            "Found msconvert, but it failed a quick check.\n"
            "Make sure this is the CLI 'msconvert.exe', not the GUI/Installer icon.\n\n"
            f"{preflight_msg}"
        )
        if log_fn:
            log_fn(err)
        raise RuntimeError(err)


    raw_path = Path(raw_path)
    out_dir = Path(dest_dir) if dest_dir is not None else (Path(out_dir) if out_dir is not None else raw_path.with_suffix(""))
    out_dir.mkdir(parents=True, exist_ok=True)

    out_name = raw_path.stem + ".mzML"
    out_mzml = out_dir / out_name
    if not overwrite and out_mzml.exists():
        i = 1
        while True:
            cand = out_dir / f"{raw_path.stem}.converted{i}.mzML"
            if not cand.exists():
                out_mzml = cand
                break
            i += 1

    cmd = [exe, str(raw_path), "--mzML", "--zlib", "--64",
           "--filter", "peakPicking true 2-",
           "--outfile", str(out_mzml.name),
           "--outdir", str(out_dir)]
    if extra_filters:
        for f in extra_filters:
            cmd += ["--filter", f]

    cp = subprocess.run(cmd, capture_output=True, text=True, check=False)

    if log_fn:
        if cp.stdout: log_fn(cp.stdout.strip())
        if cp.stderr: log_fn(cp.stderr.strip())

    if cp.returncode != 0 or not out_mzml.exists():
        raise RuntimeError(f"msconvert failed.\nSTDOUT:\n{cp.stdout}\nSTDERR:\n{cp.stderr}")
    return out_mzml
