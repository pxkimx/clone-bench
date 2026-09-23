# Clone Bench Web

The same Clone Bench interface and the same Python core (`py/core*.py`, copied from `server/`) running **inside the
browser tab** with Pyodide 0.28.3, which ships Biopython 1.85 and numpy. Nothing is uploaded: files are read by the
page and handed to Python in the tab.

    index.html    the desktop app's web/index.html with one flag (window.CB_WEB) that routes api() to Pyodide
    py/           server/core.py, core_seq.py, core_primers.py, core_sanger.py, core_protein.py
    examples/     the bundled example files

Rebuild after changing the app: `.venv/bin/python make_webapp.py` from `source/`.

## What is absent, and why
The **assistant**, the **PDF report** and **saved history** need the local server (an API key kept on your Mac, a PDF
library, a workspace folder), so this build hides them and says so on the home page. Analyses last until the tab is
closed. Every tool itself — Construct, Primers, Sanger, Protein, Sequence tools — runs here with the same code.

Biopython 1.85 differs from 1.88 in amino-acid composition (`amino_acids_percent` vs the deprecated
`get_amino_acids_percent()`, which returns fractions); `compat_aa_percent` in `py/core_protein.py` handles both, so the
numbers are the same.

## Run locally

    python3 -m http.server 8778 --directory webapp      # then open http://localhost:8778

Pyodide (~20 MB) is fetched from jsdelivr the first time an analysis runs, then cached by the browser.

## Deploy to Cloudflare Pages
Static files only — upload this folder (Workers & Pages → Create → Pages → Upload assets), or
`wrangler pages deploy . --project-name clone-bench-web`. No build step and no special headers are needed.
