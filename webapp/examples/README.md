# Example files

| File | What it is | Origin |
|---|---|---|
| `pBAD30.gb` | pBAD30 expression vector, GenBank, 4,923 bp circular (araC, pBAD promoter, bla). Its `AP(R)` CDS is annotated 3 bp past its stop codon — the Construct tool flags it. | Biopython `Tests/GenBank/pBAD30.gb` |
| `addgene-plasmid-39296-sequence-49545.gbk` | pFA6a-kanMX6 (Addgene 39296), GenBank, 3,938 bp circular, with features that cross the origin. | Biopython `Tests/GenBank/` |
| `pFA-KanMX4.dna` | pFA-KanMX4, binary SnapGene file, 3,941 bp circular. | Biopython `Tests/SnapGene/pFA-KanMX4.dna` |
| `3730.ab1` | A Sanger read from an ABI 3730 (a pCAG vector read through mCherry into a GFP fusion), with quality values. | Biopython `Tests/Abi/3730.ab1` |
| `310.ab1` | A Sanger read from an older ABI 310, **without** quality values — shows how the app reports an unscored read. | Biopython `Tests/Abi/310.ab1` |
| `demo_reference.fasta`, `demo_reference.gb` | **Made for the demo, not a real clone**: the Mott-trimmed, high-quality basecalls of `3730.ab1` with two planted edits (one substitution at 185, one extra base at 426 so the read shows a 1-bp deletion). The `.gb` also annotates the ORF the read contains (mCherry, a linker and the start of GFP) as a CDS so consequences can be shown. Rebuild with `tests/make_demo_reference.py`. | derived from `3730.ab1` |

The first five files are copied unchanged from the test suite of Biopython (https://github.com/biopython/biopython,
`Tests/`) and are redistributed under the Biopython licence (the Biopython License Agreement / BSD 3-Clause License
that covers the Biopython distribution, test files included). The plasmid sequences themselves originate from their
depositors (for example Addgene plasmid 39296) and are used here only as examples.
