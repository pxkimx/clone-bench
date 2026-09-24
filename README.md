# Clone Bench

A local app for bench cloning work: open a construct and check it, plan a digest on a virtual gel, check primers
under the buffer you actually use, and read Sanger traces against the reference — with a plain-language reading of
every figure, a PDF report, and a built-in Claude assistant. Fully offline, built on Biopython 1.88. Sister app of
RNAseq Bench and MassSpec Bench.

## Start

**Mac app** — unzip, run `xattr -cr` on the app once (see README-FIRST.txt), double-click. The first launch installs
its packages into `~/Library/Application Support/CloneBench` (about a minute); after that it needs no internet. It
serves http://localhost:8768, so it runs alongside RNAseq Bench (8765) and MassSpec Bench (8766).

**From source**
```bash
./start.sh                   # creates .venv on first run, then opens http://localhost:8768
./dev.sh                     # the same with auto-reload, for development
python -m server.selftest    # every tool on the bundled examples, offline, ~10 s
python -m pytest tests -q    # numbers checked against known values
```
The assistant needs an Anthropic API key: Settings & API key in the sidebar (stored in `~/.clone-bench/config.json`;
if none is saved there, RNAseq Bench's key is used, then `ANTHROPIC_API_KEY`).

## Tools

**Construct** — SnapGene `.dna`, GenBank, FASTA, EMBL, `.xdna`, `.gck` or a pasted sequence. Tiles (length,
topology, GC, features, single cutters); an interactive circular or linear map (features coloured by type, single
cutters around the outside, hover and click, ORF track); a features table in which every CDS is translated and
compared with its `/translation` (mismatch, internal stop, missing start/stop, length not a multiple of 3 are
flagged); restriction analysis with every commercially available enzyme (single, double, cuts-once-outside a feature
via `Analysis.do_not_cut`, only inside via `only_between`, flanking, non-cutters); a virtual digest on an agarose gel
(GeneRuler 1 kb Plus or NEB 1 kb Plus, 0.7/1/2 %, co-migration, faint bands, supercoiled uncut plasmid, comparison
with the construct minus a feature, suggested diagnostic digests); a six-frame ORF scan in any NCBI genetic code,
across the origin; exports (GenBank, FASTA, SVG/PNG figures, CSV tables, PDF report).

**Primers** — Tm by nearest-neighbour thermodynamics (`Tm_NN`, DNA_NN3, Owczarzy 2008 salt correction with Mg²⁺
and dNTPs, DMSO correction), by GC content and by the Wallace rule, side by side; MW, 3′ clamp, runs and repeats;
Taq and high-fidelity buffer presets (the latter says plainly that NEB/Thermo calculators use proprietary
parameters); a transparent complementarity heuristic for self-dimers, hairpins and cross-dimers (no ΔG — it says so);
pair analysis with an annealing suggestion; a virtual PCR on a loaded construct or pasted template (≤ 2 mismatches, a
perfect 3′ end, cloning tails, products across the origin, off-target sites) drawn on the map.

**Design primers** — for a feature or a region of a construct (or a pasted template), under the Primers tool's
buffer. *PCR*: candidates within 4 °C of the target Tm on each side, scored primer3-style (penalties from the ideal,
weights printed with the result: Tm, GC, 3′ clamp, runs, dinucleotide repeats, 3′ self-complementarity, hairpins,
cross-dimers, other binding sites on the whole template), top pairs with real alternatives. *Cloning*: primers anchored
on the insert ends and lengthened to the target Tm, with restriction-site tails (and an optional Kozak); a
start-codon site such as NcoI or NdeI is merged with the gene's ATG, and any change it forces in codon 2 is reported;
an enzyme that cuts inside the insert is flagged. *Sequencing*: a primer walk (one or both strands) with the stretch
each read should cover. Any design opens in the Primers tool for its virtual PCR; a map feature has a *Design
primers* button.

**Sanger** — `.ab1` traces with their four channels, basecalls and quality bars (zoom, scroll, jump to difference);
Mott trimming (the `abi-trim` algorithm); Q20/Q30 statistics; local alignment to the reference on both strands (as
reference + reference for circular ones); every substitution, insertion and deletion with its position, quality,
feature, codon change and consequence (synonymous, missense in 3- and 1-letter form, nonsense, frameshift) and a
confidence call (confident / check the trace / probably a basecall error); mixed peaks reported separately; a
coverage track per read and combined, answering “is my clone fully sequence-verified?” for any feature. The demo is
labelled for what it is: the reference is the read's own basecalls with two planted edits.

**Protein** — ProtParam properties (MW, pI, charge at pH 7.4 and a charge–pH curve, ε280 reduced and with cystines,
GRAVY, instability index with its caveat, aromaticity, composition), an A280 → mg/mL and µM converter (flags a
protein without Trp/Tyr), Kyte–Doolittle hydropathy (window 9 / 19 with the TM threshold), tags (His6, FLAG, HA,
Myc, Strep-tag II, V5, GST/MBP starts) and protease sites (TEV, HRV 3C, thrombin, enterokinase) with their scars.

**Sequence tools** — reverse complement, six-frame translation in any NCBI code, IUPAC motif search on both strands
(across the origin for circular sequences), SEGUID checksums with a rotation- and strand-independent circular form.

## Browser version
`webapp/` is the same UI and the same Python core (`server/core*.py`) running inside Pyodide 0.28.3 (Biopython
1.85), built by `make_webapp.py`. The assistant and the PDF report need the local server and are absent there; the page says so. Analyses are saved
in the browser (IndexedDB) and can be downloaded to a `.json` file and opened again. Serve it with `python -m http.server 8778 --directory webapp`. See `webapp/README.md`.

## Layout
```
server/core.py           api(name, args): the one entry point    server/app.py      HTTP API (+ workspace)
server/core_seq.py       reading, CDS checks, restriction,       server/common.py   workspace, settings
                         digest/gel, ORFs, exports, seq tools    server/report.py   PDF (reportlab)
server/core_primers.py   Tm, heuristics, virtual PCR             server/agent.py    Claude assistant
server/core_sanger.py    traces, trimming, alignment, calls      server/selftest.py offline self-test
server/core_protein.py   ProtParam, A280, tags                   web/index.html     the whole UI
server/core_design.py    primer design (PCR, cloning, sequencing)
```
Workspace: `~/Library/Application Support/CloneBench/work/<id>/` (inputs, `request.json`, `result.json`,
`report.pdf`). Set `CB_HOME` to move it.

## Caveats
The dimer/hairpin check finds complementarity, not stability — use primer3 or a vendor tool for ΔG. Vendor Tm
calculators differ by a few °C. The gel is an illustration, not a calibrated model. Methylation (Dam/Dcm/CpG)
blocking of restriction sites is not modelled. A Sanger difference below Q20 is reported as probably not real.
Everything runs on your machine; only the assistant talks to Anthropic's API, and only with what you type plus the
tool results it requests.
