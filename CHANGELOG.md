# Changelog

## 0.2.0 — 2026-09-23
- **New tool: Design primers.** Clone Bench could check primers but not design them. Pick a construct and a feature
  (or a region, or paste a template) and one of three modes. All of them use the Primers tool's buffer and the same
  Tm, clamp, dimer, hairpin and off-target checks, so a designed primer is judged exactly as a checked one would be.
  - **PCR:** the best pairs amplifying the target. They are scored like primer3 (penalties from the ideal, with the
    weights printed under the table) and checked for cross-dimers and for other binding sites on the whole
    template. No primer is used in more than two listed pairs, so the list offers real alternatives.
  - **Cloning:** primers anchored on the insert ends, lengthened to the target Tm, with restriction-site tails and
    an optional Kozak sequence.
    - A start-codon site (NcoI CCATGG, NdeI CATATG) is merged with the gene's ATG instead of being put in front
      of it. Putting it in front leaves an out-of-frame ATG ahead of the gene in a pET-type vector.
    - When the site forces a change in codon 2, the tool says so (for example "K2E, an amino-acid change; NdeI
      needs none").
    - An enzyme that cuts inside the insert is flagged.
  - **Sequencing:** a primer walk on one or both strands, with the stretch each read should cover drawn on the map.
- Any design opens in the Primers tool (virtual PCR, dimers) with one click. Clicking a feature on a construct map
  now offers **Design primers**. The assistant has a `design_primers` tool.
- Tested the way the other tools are: about 950 randomised designs on generated plasmids, checked against the right
  answer. Each PCR pair was verified to sit on the template at the stated positions, give the stated product in a
  virtual PCR, cover the whole target and carry the exact Tm. Each cloning primer was checked to match the insert
  ends, keep the first ATG in frame, and warn exactly when an enzyme cuts the insert. Each sequencing walk was
  checked to cover its target.
- **Fixed in the browser version.** The page kept its own list of Python modules, separate from the build script's,
  so a new module shipped but was never loaded. The build now writes the list into the page.

## 0.1.1 — 2026-09-23
Fixes found by randomised testing, which ran each tool thousands of times on inputs with a known right answer:
random plasmids with planted genes, primers cut from known positions, and a real Sanger read against references
with edits planted in them. The answers were checked against Biopython directly.
- **A primer named R, or KanR, gained extra "bases".** "R ACGT…" and "KanR TTAG…" are how people name primers, but
  every letter of R and of KanR is also an IUPAC code, so the name was read as the primer's first bases. The primer
  came out 1–4 nt long, the virtual PCR product 1–4 bp long, and one mismatch was reported that isn't there. The
  bundled self-test used exactly this pair and still passed, because the extra bases looked like a cloning tail. A
  first word in mixed case, or of three letters or fewer ahead of a real sequence, is now a name.
- **No Sanger PDF report could be written.** Its methods text says "Q < 20". The PDF writer read the "<" as the start
  of a tag and deleted text up to the next ">", which took part of a formatting tag with it, so ReportLab stopped.
  Bare "<" is now escaped. A paragraph that still can't be formatted goes in as plain text rather than failing the
  report.
- **Five enzymes that never cut a plasmid prep were listed.** AbaSI, FspEI, LpnPI, MspJI and SgeI cut only
  methylated or hydroxymethylated DNA. Searched as plain patterns, they "cut" at almost every C and filled the enzyme
  table with nonsense. They are left out now, and the page says so.
- **A FASTA with a pasted symbol crashed.** A curly quote or Greek letter copied from Word or a PDF gave an internal
  error. You now get a sentence naming the odd characters.
- **A buffer with more dNTPs than Mg²⁺ and no Na⁺/K⁺/Tris crashed.** dNTPs bind Mg²⁺, so no free cations were left
  and the salt correction divided by zero. You now get an explanation. Whenever Mg²⁺ ≤ total dNTPs there is also a
  warning: no Mg²⁺ is free for the polymerase, and the usual cause is entering dNTPs per nucleotide rather than total.
- **Reports download as `date_toolname_report`**, for example `2026-09-23_CloneBench-Sanger_report.pdf`. The date is
  the day the analysis ran.

## 0.1.0 — 2026-09-23
First release.
- **Construct**: SnapGene/GenBank/FASTA/EMBL/xdna/gck or pasted sequence; circular and linear interactive maps;
  features table with CDS translation checks against /translation; restriction analysis with every commercially
  available enzyme and feature-aware filters (`do_not_cut`, `only_between`, flanking, non-cutters); virtual digest
  and agarose gel with two real ladders, co-migration, supercoiled uncut plasmid, "construct minus feature"
  comparison and suggested diagnostic digests; six-frame ORFs in any NCBI code across the origin; GenBank/FASTA/
  SVG/PNG/CSV exports.
- **Primers**: Tm_NN (DNA_NN3, Owczarzy 2008 salt correction, DMSO) beside Tm_GC and Wallace; buffer presets kept for
  the session; clamp, runs, repeats; complementarity heuristic for dimers/hairpins (labelled as such); pair analysis;
  virtual PCR with cloning tails, off-targets and products across the origin.
- **Sanger**: trace viewer, Mott trimming, both-strand local alignment, every difference with quality, feature,
  codon consequence and a confidence call, mixed peaks, coverage/verification per feature; honest planted-edit demo.
- **Protein**: ProtParam properties, charge–pH curve, A280 → concentration, hydropathy, composition, tags and protease
  scars; Biopython 1.85/1.88 composition shim.
- **Sequence tools**: reverse complement, translation, IUPAC motifs, SEGUID with a circular (rotation/strand-free) form.
- Workspace history, PDF report, in-app Claude assistant with ten tools, offline self-test, pytest suite.
- **Browser build** (`webapp/`): the same UI and core in Pyodide 0.28.3 / Biopython 1.85, via `make_webapp.py`.
