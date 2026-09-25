#!/usr/bin/env python3
"""Build the browser-only version into webapp/ from the same UI and the same Python core.

    .venv/bin/python make_webapp.py

- web/index.html → webapp/index.html, with one flag injected (window.CB_WEB = true) that makes the page's single
  api(name, args) function run the Python core inside Pyodide instead of POSTing to /api/run/<name>
- server/core*.py → webapp/py/ (checked here to import nothing but Biopython, numpy and the standard library)
- the example files → webapp/examples/

The assistant needs the local server, so the page hides it and says so. The desktop PDF report (reportlab) is not
available either, but the browser build has its own "Download report" button that writes a self-contained HTML
report instead — see webReport() in web/index.html.
Serve locally with:  .venv/bin/python -m http.server 8778 --directory webapp
"""
from __future__ import annotations

import ast
import shutil
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "webapp"
CORE = ["core.py", "core_seq.py", "core_primers.py", "core_sanger.py", "core_protein.py", "core_design.py"]
EXAMPLES = ["pBAD30.gb", "addgene-plasmid-39296-sequence-49545.gbk", "pFA-KanMX4.dna", "3730.ab1", "310.ab1",
            "demo_reference.gb", "demo_reference.fasta", "README.md"]
ALLOWED_TOP = {"Bio", "numpy", "__future__"}


def check_imports(path: Path):
    """The browser build can only run code that needs Biopython, numpy and the standard library."""
    tree = ast.parse(path.read_text())
    stdlib = set(sys.stdlib_module_names)
    bad = []
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names = [node.module.split(".")[0]]
        for n in names:
            if n not in ALLOWED_TOP and n not in stdlib:
                bad.append(n)
    if bad:
        raise SystemExit(f"{path.name} imports {sorted(set(bad))}: not allowed in the core (Biopython, numpy, stdlib only).")


def main():
    (OUT / "py").mkdir(parents=True, exist_ok=True)
    (OUT / "examples").mkdir(parents=True, exist_ok=True)
    for f in CORE:
        src = ROOT / "server" / f
        check_imports(src)
        shutil.copy2(src, OUT / "py" / f)
    for f in EXAMPLES:
        shutil.copy2(ROOT / "examples" / f, OUT / "examples" / f)
    html = (ROOT / "web" / "index.html").read_text()
    flag = ("<script>window.CB_WEB = true;  /* browser build: api() runs the Python core in Pyodide */\n"
            f"window.CB_CORE = {json.dumps(CORE)};</script>\n")
    marker = "<meta charset=\"utf-8\">\n"
    assert marker in html, "web/index.html changed its <head>; update make_webapp.py"
    html = html.replace(marker, marker + flag, 1)
    html = html.replace("<title>Clone Bench</title>", "<title>Clone Bench Web</title>", 1)
    (OUT / "index.html").write_text(html)
    print(f"webapp/ built: index.html, {len(CORE)} core modules, {len(EXAMPLES)} example files")


if __name__ == "__main__":
    main()
