# OrthoFinder -> Ortholog Extractor

This repository contains two small tools for extracting ortholog candidates from OrthoFinder output:
- `extract_from_orthogroups.py` — command-line tool (CLI)
- `extract_from_orthogroups_gui.py` — Tkinter graphical tool (GUI)

Purpose
-------
These tools read the `Orthogroups.tsv` output from OrthoFinder and a CSV of query IDs, and extract orthogroup members (per species) for each query. They support matching either full protein IDs or gene-prefixes (e.g., `AAEL027978` matching `AAEL027978-PA`) and can output either a single wide table across species or per-target species outputs. They also write a "not found" file listing queries that did not match any orthogroup.

Credit — OrthoFinder
---------------------
These scripts operate on OrthoFinder output (they DO NOT run OrthoFinder themselves). OrthoFinder is an independent tool — please cite and credit it when you use its output:

- OrthoFinder (Emms & Kelly). Project: https://github.com/davidemms/OrthoFinder  
  Paper / citation: Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology.

Prerequisites
-------------
- Python 3.8+ (the scripts use only the Python standard library).
- For the GUI: Tk / Tcl (Tkinter) must be available in your Python environment. If Tkinter is not present, use the CLI.
  - Debian/Ubuntu: `sudo apt install python3-tk`
  - Conda: `conda install -c conda-forge tk`
  - macOS/Windows: usually included with standard Python installers; use your distribution's instructions if missing.
- Optional: `--id-map` supports a two-column TSV mapping protein_id -> gene_name; this file is optional and no extra Python packages are required.

Generate OrthoFinder output
---------------------------
1. Install OrthoFinder and run it on your set of proteomes as described in its documentation:
   - OrthoFinder repository and docs: https://github.com/davidemms/OrthoFinder
2. After OrthoFinder finishes, locate the `Orthogroups.tsv` file (usually in the `Results_<date>/` directory produced by OrthoFinder).

Typical OrthoFinder run (very short example)
- Prepare FASTA files (one per species).
- Run OrthoFinder (example):
  - `orthofinder -f /path/to/proteomes/ -t 8 -a 8`
- After run, note the path to `Orthogroups.tsv` in the results folder.

What these scripts do
---------------------
- Parse `Orthogroups.tsv` (first column is Orthogroup ID, subsequent columns are species).
- Build an index mapping protein IDs to orthogroups and a gene‑prefix index so gene-prefix matching works (handles common separators like `-`, `|`, `.`, `_`, `:`).
- Read your CSV of queries, extract IDs from a selected column (or a newline list in older CLI usage).
- For each query:
  - Find matching proteins in the specified source species (exact or gene-prefix).
  - For matched orthogroup(s), output orthogroup members for the chosen target species (or produce a wide table across all species, excluding the source species when comparing to `All`).
  - Write a separate "not found" file listing queries that had no matches.

Files written
-------------
- If target is a specific species:
  - `<out_prefix>_<target>_target.tsv` — columns: query_id, orthogroup, `<target>_members`
  - `<out_prefix>_<target>_not_found.tsv` — queries with no orthogroup match (query_id, reason)
- If target == `All`:
  - `<out_prefix>_all_wide.tsv` — columns: query_id, orthogroup, one column per species (source species excluded)
  - `<out_prefix>_all_not_found.tsv` — queries with no orthogroup match

GUI usage (Tkinter)
-------------------
1. Go to where the repository is located e.g :
   - `cd orthogroup_extractor`
2. Run:
   - `python extract_from_orthogroups_gui.py`
3. In the GUI:
   - Paste full path to `Orthogroups.tsv` and click "Load species".
   - Paste full path to your CSV file (or type it) and click "Load CSV columns".
   - Select the CSV column containing query IDs.
   - Choose the source species (where your query IDs come from).
   - Choose a target species or pick `All`.
   - (Optional) Toggle "Match gene-prefix" — the GUI will auto-enable it if the CSV column looks like gene IDs.
   - Set output folder and prefix, then click "Run".
4. Inspect the written files in the output folder and the GUI log / popup summary.

Notes and troubleshooting
-------------------------
- Exact species names: species names come from `Orthogroups.tsv` header — use the exact string (case and punctuation matter).
- Gene-prefix matching: many proteomes encode isoform suffixes (e.g., `AAEL027978-PA`); the tools strip common separators and build a gene-prefix index so queries like `AAEL027978` match all isoforms.
- If GUI fails with a tkinter import error, install system Tk or use the CLI (the CLI does not require Tkinter).
- The scripts are designed to avoid external pip/conda dependencies; if you removed the optional `pyperclip` import block, no third-party package is required.

Repository layout suggestion
---------------------------
```
.
├── README.md
├── extract_from_orthogroups.py
├── extract_from_orthogroups_gui.py
└── sample_data/               # optional: small example Orthogroups.tsv + queries.csv
```

License
-------
MIT License

Copyright (c) 2025 bertiebertolo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Acknowledgements
----------------
- OrthoFinder (Emms & Kelly) for orthogroup inference; link: https://github.com/davidemms/OrthoFinder  
- This repository contains utilities to extract orthogroup information from OrthoFinder output for downstream ortholog discovery.
