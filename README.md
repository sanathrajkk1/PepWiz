# PepWiz 🧪  
**A GUI-based peptide MS/MS processing and visualization tool**

PepWiz helps you process peptide mass spectrometry (MS/MS) data with ease — no coding required.  
Provide a peptide sequence and your RAW or mzML/mzXML file, and PepWiz will automatically convert, filter, match, and visualize fragments, generating clean, editable figures for papers or presentations.

---

## ✨ Features

- 🧬 **Automatic b/y fragment matching** — theoretical vs. observed peaks  
- ⚙️ **ProteoWizard integration** — optional RAW → mzML conversion via `msconvert`  
- 🎨 **Publication-ready visuals** — editable SVG fragment map & annotated spectrum  
- 🪄 **Intuitive GUI** — one-click processing, no scripts needed  
- 💾 **Exports**
  - `*.out` — text summary  
  - `*.fragments.svg` — sequence coverage map  
  - `*.spectrum.svg` — annotated spectrum  

---

## 🧩 Installation Guide

PepWiz runs on **Windows with Python 3.10+**.  
ProteoWizard’s `msconvert` is optional but required if you plan to open vendor `.RAW` files.

### 🟢 Option 1 — Install directly from GitHub (recommended)

```bash
pip install git+https://github.com/sanathrajkk1/PepWiz.git
```
Then launch:
```bash
python -m pepwiz.gui
```

### 🟢 Option 2 — Local install (manual or offline setup)
```bash
git clone https://github.com/sanathrajkk1/PepWiz.git

cd pepwizard

pip install .
```
#(For development)
```bash
python -m venv .venv

.\.venv\Scripts\activate

pip install -e .[io]

python -m pepwiz.gui
```
---

### ⚙️ ProteoWizard Setup (for .RAW files)

If you only work with .mzML or .mzXML, skip this section.

- Step 1 — Install ProteoWizard
	Download from: https://proteowizard.sourceforge.io/download.html

- Step 2 — Find the real msconvert.exe
Don’t use the GUI shortcut (MSConvertGUI_Icon.exe).

You need the command-line binary.

```cmd
	where /r C:\ msconvert.exe
```
```powershell
	Get-ChildItem "$env:LOCALAPPDATA","$env:ProgramFiles","$env:ProgramFiles(x86)" -Recurse -Filter msconvert.exe -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty FullName
```
- Step 3 — Register the path
```cmd
	setx PEPWIZ_MS_CONVERT "C:\Program Files\ProteoWizard\msconvert.exe"
```
(Replace with your actual path.)
Restart your terminal — PepWiz will detect it automatically.

---

### 🚀 Using PepWiz

Launch:
```cmd
	python -m pepwiz.gui
```
Workflow:

1. Open MS File → choose .raw, .mzML, or .mzXML

2. Enter your peptide sequence

3. Adjust:

        PPM tolerance
        Charge states
        Optional RT window

4. Click Run and monitor logs

Outputs:

    filename.out — matched fragments
    filename.fragments.svg — sequence coverage
    filename.spectrum.svg — annotated spectrum

---

### 🧪 Developer & Contributor Setup
```cmd
git clone https://github.com/sanathrajkk1/PepWiz.git
cd pepwiz
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .[io]
pytest -q
```
---

### 🧠 Troubleshooting

"PepWiz can’t find msconvert"

1. Install ProteoWizard

2. Find binary
```cmd
		where /r C:\ msconvert.exe
```
3. Register path
```cmd
		setx PEPWIZ_MS_CONVERT "C:\path\to\msconvert.exe"
```
4. Restart and rerun PepWiz.

"msconvert failed or DLL error"

- Reinstall full ProteoWizard CLI (not GUI-only).
	
- Ensure you reference msconvert.exe, not MSConvertGUI_Icon.exe.

"GUI not opening"

- Python ≥ 3.10

- tkinter installed (default on Windows)

- Run via terminal:
```cmd
		python -m pepwiz.gui
```
---

### 🧬 Example Run

1. Convert or download a small .mzML file 

2. Run:
```cmd
		python -m pepwiz.gui
```		

3. Parameters:

        Sequence: "PEPTIDE"
        PPM: 10
        Charge: 2

4. Click Run

5. Outputs will appear in the same folder.
	
---
	
### 👩‍🔬 Contributors & Contact

- Sanath Raj Kavuthian Kandy — Concept, design, GUI development, optimization, validation

- Jonathan Chekan — Concept, validation

📧 sanathrajkk1@gmail.com
Feedback and pull requests welcome!

---

### 📚 Citation

    

---

### ⚖️ License

This project is licensed under the MIT License.
See LICENSE for full details.

---