# Changelog

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
