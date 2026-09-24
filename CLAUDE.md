# Clone Bench — notes for Claude Code

Local app for bench cloning: construct maps and checks, restriction/digest/gel, primers, Sanger verification,
protein properties, sequence tools. Sister of RNAseq Bench and MassSpec Bench (same look: FastAPI + one-file
vanilla-JS UI + in-app Claude assistant + PDF report). Owner: Paul (bench biologist). Every figure has a **how**
(what it shows, how computed) and a **yours** (what it says about THIS input, generated from the numbers).

## Architecture — keep it
- All computation lives in `server/core*.py`, which import **only Biopython, numpy and the standard library**
  (`make_webapp.py` enforces this with an import check). Public entry point: `core.api(name, args) -> dict`, JSON in
  and out. `UserFacingError` = a message for the user, shown without a traceback.
- The UI calls exactly one function, `api(name, args)` in `web/index.html`: locally a POST to `/api/run/<name>`; in the
  browser build (`window.CB_WEB`) the same core inside Pyodide. Never add a UI feature that bypasses `api()` for
  computation — it would silently break the web build.
- `server/app.py` saves the whole-tool analyses (construct, primers, design, sanger, protein, seqtools) as workspace
  items under `CB_HOME/work/<id>/` (inputs + request.json + result.json). Follow-ups (digest, enzymes, orfs,
  suggest_digest, export, translate_feature) are computed on demand from the item's `record` and not saved.
- Result format (`core.Result`): tiles, flags (ok/warn/info), sections of items — `fig` (drawn in the browser by
  `DRAW[draw]` from `data`), `table`, `widget` (`WIDGET[widget]`, interactive, calls `api`) — methods, versions.
  Constructs, primer templates and Sanger references carry `record` (sequence + features JSON; parts in biological
  order, 0-based half-open, `start0/end0/wraps` for the outer span).
- Figures are SVG/canvas drawn client-side with CSS variables; `svgString()` swaps `var(--x)` for the light palette
  (`LIGHT`) on export, so downloads and the PDF are always light. The PDF gets PNG snapshots of figures (`card.__snap`)
  and current widget tables (`card.__tables`) POSTed to `/api/items/<id>/report`.
- Assistant (`server/agent.py`): same SSE protocol as RNAseq Bench; tools call `core.api` on saved items.

## Conventions / gotchas
- Coordinates shown to users are 1-based inclusive. Restriction cut positions are Biopython's convention: 1-based
  first base after the top-strand cut (`search()`); internally `pos0 = pos - 1` is the cut coordinate.
- `Bio.Restriction` keeps search state on the enzyme classes: every search goes through `core_seq._RLOCK`.
- Circular: restriction via `linear=False`; ORFs and alignment scan sequence + sequence; features across the origin
  are compound locations whose parts end at L and start at 0 (`feature_span`).
- Biopython 1.85 (Pyodide) vs 1.88: `compat_aa_percent` in core_protein.py. Pyodide's 1.85 already has the
  `amino_acids_percent` property; its `get_amino_acids_percent()` returns fractions. Test both in tests/.
- Honesty: dimers/hairpins are a heuristic (say so), vendor Tm differ, gel is illustrative, Sanger diffs at Q < 20 are
  "probably a basecall error", the demo reference is labelled as planted. Don't add a number without its caveat.

- Primer design (`core_design.py`) reuses core_primers' `tm_nn`, `buffer_from`, `comp_scan`, `hairpin_scan` and
  `binding_sites` — a designed primer must be judged exactly like a checked one. Circular targets work in a rotated
  `Frame` (`.orig()` maps back). Start-codon sites (NcoI, NdeI) are merged with the ATG; never prepend a whole site.
- Browser build saving: `savedPut/savedLoad/savedDownload/savedOpen/savedClear` in web/index.html keep whole-tool
  results (`SAVE_KINDS`) in IndexedDB `clone-bench-web` and a `{app:"Clone Bench",format:1,items}` .json file. The
  stored copy is taken by `put()` inside `api()`, before the UI mutates the result. Keep new result kinds in step
  with `SAVE_KINDS` and the list in `api()`.
- The browser build loads the core modules listed in `make_webapp.CORE`, injected into the page as `window.CB_CORE`.
  Add a new `core_*.py` there (only) — a hand-kept second list in index.html once shipped a module without loading it.

## Run / build / test
- `./dev.sh` (port 8768, auto-reload) · `python -m server.selftest` · `python -m pytest tests -q`
- `python make_webapp.py` → `webapp/` (serve: `python -m http.server 8778 --directory webapp`)
- `python tests/make_demo_reference.py` rebuilds `examples/demo_reference.{fasta,gb}`
- `./build_mac.sh` → `../Clone Bench.app` + zip; `macos/make_icon.py` draws the icon from the logo geometry.
- Bump `VERSION` + `CHANGELOG.md` for every change that ships; the launcher uses VERSION to replace a running older server.
