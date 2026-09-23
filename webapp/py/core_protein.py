"""Protein: ProtParam properties, A280 → concentration, hydropathy, composition, tags and protease sites.
Biopython + stdlib only.

Amino-acid composition differs between Biopython releases: the cached property `amino_acids_percent` returns
0–100, while the older method `get_amino_acids_percent()` returns fractions (and in releases before the property
existed, `amino_acids_percent` was a plain attribute that stayed None until that method had been called). Checked in
the browser build: Pyodide's Biopython 1.85 already has the property (percent) and a deprecated
`get_amino_acids_percent()` (fractions). `compat_aa_percent` prefers the property, falls back to the method and
rescales fractions to percent, so every version gives the same numbers; tests/test_core.py exercises both paths.
"""
from __future__ import annotations

import functools
import re

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

from .core import Result, UserFacingError, esc, pl
from .core_seq import extract_feature, find_feature, normalize_dna

STD = "ACDEFGHIKLMNPQRSTVWY"

TAGS = [  # name, regex, note
    ("His-tag", r"H{6,}", "binds Ni-NTA / Co resin"),
    ("FLAG", r"DYKDDDDK", "anti-FLAG M2; its last five residues are also an enterokinase site"),
    ("HA", r"YPYDVPDYA", "anti-HA"),
    ("Myc", r"EQKLISEEDL", "anti-Myc 9E10"),
    ("Strep-tag II", r"WSHPQFEK", "Strep-Tactin"),
    ("V5", r"GKPIPNPLLGLDST", "anti-V5"),
    ("GST (start)", r"MSPILGYWKIKGLVQPTRLL", "the first 20 residues of S. japonicum GST, as in pGEX (heuristic: matches the start only)"),
    ("MBP (start)", r"KIEEGKLVIWINGDKGYNG", "the N-terminal stretch of E. coli MalE, as in pMAL (heuristic: matches the start only)"),
]
PROTEASES = [  # name, regex, cut offset within the match (residues before the cut), scar note
    ("TEV", r"ENLYFQ[GS]", 6, "cuts ENLYFQ↓G/S: the downstream protein keeps one extra N-terminal G or S"),
    ("HRV 3C (PreScission)", r"LEVLFQGP", 6, "cuts LEVLFQ↓GP: the downstream protein keeps GP at its N-terminus"),
    ("Thrombin", r"LVPRGS", 4, "cuts LVPR↓GS: the downstream protein keeps GS at its N-terminus"),
    ("Enterokinase", r"DDDDK(?!P)", 5, "cuts DDDDK↓: no residues left on the downstream protein (not before a proline)"),
]


def compat_aa_percent(pa) -> dict:
    """Percent (0–100) of each of the 20 amino acids, on Biopython 1.85 and on 1.87+."""
    attr = getattr(type(pa), "amino_acids_percent", None)
    if isinstance(attr, (property, functools.cached_property)):
        v = pa.amino_acids_percent                 # 1.87+: a (cached) property
    elif callable(getattr(pa, "get_amino_acids_percent", None)):
        v = pa.get_amino_acids_percent()           # 1.85: a method returning fractions
    else:                                          # neither: count it ourselves
        seq = str(pa.sequence)
        v = {a: seq.count(a) / max(1, len(seq)) for a in STD}
    v = dict(v)
    total = sum(v.values())
    scale = 100.0 if total <= 1.5 else 1.0         # fractions → percent; percent stays percent
    return {a: float(v.get(a, 0.0)) * scale for a in STD}


def clean_protein(text: str) -> tuple[str, list]:
    notes = []
    s = re.sub(r"[\s\d]+", "", text or "").upper()
    if s.startswith(">"):
        raise UserFacingError("Paste the sequence without its FASTA header line, or drop the file.")
    if not s:
        raise UserFacingError("The protein sequence is empty.")
    if s.endswith("*"):
        s = s[:-1]
    if "*" in s:
        k = s.index("*")
        notes.append(("warn", f"A stop (*) at residue {k + 1}; only the {k} residues before it were analysed."))
        s = s[:k]
    odd = {c: s.count(c) for c in set(s) - set(STD)}
    if odd:
        bad = set(odd) - set("XBZJUO-.")
        if bad:
            raise UserFacingError(f"The sequence contains characters that are not amino acids: {', '.join(sorted(bad))}.")
        s = "".join(c for c in s if c in STD)
        notes.append(("warn", f"{sum(odd.values())} non-standard residue{'s' if sum(odd.values()) > 1 else ''} "
                              f"({', '.join(f'{k}×{v}' for k, v in sorted(odd.items()))}) were left out: ProtParam knows only the "
                              "20 standard amino acids" + (" (U is selenocysteine)" if "U" in odd else "") + "."))
    if len(s) < 5:
        raise UserFacingError("The protein is shorter than 5 residues — too short for these calculations.")
    return s, notes


def looks_like_dna(s: str) -> bool:
    t = re.sub(r"[\s\d]+", "", s or "").upper()
    return bool(t) and sum(t.count(c) for c in "ACGTUN") / len(t) > 0.9


def translate_dna(text: str, table: int) -> tuple[str, list]:
    dna, notes = normalize_dna(text, "The DNA")
    try:
        ct = CodonTable.unambiguous_dna_by_id[int(table)]
    except (KeyError, ValueError):
        raise UserFacingError(f"NCBI genetic code {table} does not exist.")
    extra = len(dna) % 3
    if extra:
        notes.append(("warn", f"The DNA is {len(dna)} bp, not a multiple of 3; the last {extra} base{'s' if extra > 1 else ''} were ignored."))
        dna = dna[:len(dna) - extra]
    p = str(Seq(dna).translate(table=int(table)))
    if dna[:3] in ct.start_codons and p and p[0] != "M":
        p = "M" + p[1:]
        notes.append(("info", f"Starts with {dna[:3]}, read as Met (a valid start in table {table})."))
    elif dna[:3] not in ct.start_codons:
        notes.append(("info", f"The DNA does not start with a start codon ({dna[:3]}); translated from base 1 anyway."))
    return p, notes


def segments(values, start_pos, thr):
    """Stretches where a sliding-window average exceeds thr: [(start, end)] as 1-based residue centres."""
    out, cur = [], None
    for i, v in enumerate(values):
        if v > thr and cur is None:
            cur = i
        if v <= thr and cur is not None:
            out.append((cur + start_pos, i - 1 + start_pos))
            cur = None
    if cur is not None:
        out.append((cur + start_pos, len(values) - 1 + start_pos))
    return out


def protein(args: dict) -> dict:
    notes = []
    mode = args.get("mode") or "auto"
    table = int(args.get("table") or 1)
    name = args.get("name") or "Protein"
    src = "pasted protein"
    if args.get("record") and args.get("feature") not in (None, ""):
        rj = args["record"]
        f = find_feature(rj, args["feature"])
        if f is None:
            raise UserFacingError(f"No feature called “{args['feature']}” in {rj.get('name')}.")
        q = f.get("qualifiers") or {}
        table = int(q.get("transl_table") or args.get("table") or 11)
        seq, notes = translate_dna(extract_feature(rj["seq"], f)[int(q.get("codon_start") or 1) - 1:], table)
        name, src = f"{f['name']} ({rj.get('name')})", f"CDS {f['name']} translated with NCBI table {table}"
    else:
        text = args.get("text") or ""
        if text.lstrip().startswith(">"):
            lines = text.strip().splitlines()
            name = lines[0][1:].strip()[:60] or name
            text = "\n".join(lines[1:])
        if mode == "dna" or (mode == "auto" and looks_like_dna(text)):
            seq, notes = translate_dna(text, table)
            src = f"DNA translated with NCBI table {table}"
            if mode == "auto":
                notes.insert(0, ("info", "The input looked like DNA, so it was translated. Choose “protein” if it really is a "
                                         "protein made of A, C, G and T."))
        else:
            seq = text
    seq, n2 = clean_protein(seq)
    notes += n2
    pa = ProteinAnalysis(seq)
    L = len(seq)
    mw = pa.molecular_weight()
    pi = pa.isoelectric_point()
    ch74 = pa.charge_at_pH(7.4)
    red, ox = pa.molar_extinction_coefficient()
    nW, nY, nC = seq.count("W"), seq.count("Y"), seq.count("C")
    gravy = pa.gravy()
    ii = pa.instability_index()
    arom = pa.aromaticity()
    comp = compat_aa_percent(pa)
    curve = [[round(ph / 10, 1), round(pa.charge_at_pH(ph / 10), 3)] for ph in range(20, 121, 2)]
    hyd = {}
    for w in (9, 19):
        if L >= w:
            vals = pa.protein_scale(kd, w, 1.0)
            hyd[w] = {"start": w // 2 + 1, "values": [round(v, 3) for v in vals]}
    tm = segments(hyd[19]["values"], hyd[19]["start"], 1.6) if 19 in hyd else []

    R = Result("protein", name)
    R.r["protein"] = {"seq": seq, "mw": mw, "ext_reduced": red, "ext_cystines": ox}
    R.tile(f"{L:,} aa", "length")
    R.tile(f"{mw / 1000:,.2f} kDa", "average mass")
    R.tile(f"{pi:.2f}", "isoelectric point")
    R.tile(f"{ch74:+.1f}", "net charge at pH 7.4")
    R.tile(f"{red:,}", "ε280 (M⁻¹cm⁻¹, reduced)")
    R.r["summary"] = {"length": L, "mw": mw, "pi": pi}
    for lvl, t in notes:
        R.flag(lvl, t)
    if nW == 0 and nY == 0:
        R.flag("warn", "<b>No Trp and no Tyr</b>: this protein barely absorbs at 280 nm, so A280 cannot measure its concentration. "
                       "Use a dye (Bradford/BCA) or A205.")
    elif nW == 0:
        R.flag("info", f"No Trp, {nY} Tyr: the 280 nm signal is weak ({red:,} M⁻¹cm⁻¹), so small contaminants weigh more.")
    else:
        R.flag("ok", f"{nW} Trp and {nY} Tyr: ε280 = {red:,} M⁻¹cm⁻¹ reduced, {ox:,} with all Cys paired; A280 is a fair way to measure it.")
    R.flag("info" if ii < 40 else "warn", f"Instability index {ii:.1f} ({'below' if ii < 40 else 'above'} 40: predicted "
                                           f"{'stable' if ii < 40 else 'unstable'}). A 1990 dipeptide heuristic (Guruprasad et al.) "
                                           "— treat it as a hint, not a measurement.")
    if tm:
        R.flag("info", f"{pl(len(tm), 'stretch', 'stretches')} with a 19-residue hydropathy average above 1.6 — candidate "
                       "transmembrane helices: " + ", ".join(f"{a}–{b}" for a, b in tm[:6]) + ".")

    # ---- properties
    R.section("props", "ProtParam", "Physico-chemical properties", f"From {esc(src)}; {L} residues.")
    abs01_red, abs01_ox = red / mw, ox / mw
    R.table("props", "Properties", ["Property", "Value", "Note"], [
        ["Average molecular weight", f"{mw:,.1f} Da", "average isotopic masses; no modifications"],
        ["Isoelectric point", f"{pi:.2f}", "Bjellqvist pK values"],
        ["Net charge at pH 7.4", f"{ch74:+.2f}", ""],
        ["ε280, reduced", f"{red:,} M⁻¹cm⁻¹", f"{nW} Trp × 5,500 + {nY} Tyr × 1,490"],
        ["ε280, cystines", f"{ox:,} M⁻¹cm⁻¹", f"+ {nC // 2} cystine × 125 (all {nC} Cys paired)"],
        ["Abs 0.1 % (1 mg/mL)", f"{abs01_red:.3f} (reduced) · {abs01_ox:.3f} (cystines)", "A280 of a 1 mg/mL solution, 1 cm"],
        ["GRAVY", f"{gravy:+.3f}", "mean Kyte–Doolittle hydropathy; > 0 is hydrophobic"],
        ["Instability index", f"{ii:.1f}", "< 40 predicted stable (Guruprasad et al. 1990 heuristic)"],
        ["Aromaticity", f"{arom:.3f}", "fraction F + W + Y (Lobry 1994)"],
    ], how="Biopython <code>Bio.SeqUtils.ProtParam.ProteinAnalysis</code>, the ExPASy ProtParam methods.",
        yours=f"{mw / 1000:,.1f} kDa, pI {pi:.2f}, net charge {ch74:+.1f} at pH 7.4, GRAVY {gravy:+.2f}.", wide=False)
    R.fig("charge", "Charge vs pH", "curve", {"points": curve, "x": "pH", "y": "net charge", "mark": [7.4, ch74], "zero": True,
                                              "pi": pi},
          how="Net charge from the Henderson–Hasselbalch equation over every ionisable group (<code>charge_at_pH</code>), pH 2–12. "
              "It crosses zero at the pI; proteins are least soluble near their pI.",
          yours=(f"Net charge {ch74:+.1f} at pH 7.4, crossing zero at pH {pi:.2f}. "
                 + ("For ion exchange at pH 7.4 it binds an anion exchanger (Q)." if ch74 < -1 else
                    "For ion exchange at pH 7.4 it binds a cation exchanger (S)." if ch74 > 1 else
                    "It is nearly neutral at pH 7.4, so ion exchange needs a pH further from the pI.")))

    # ---- A280
    R.section("a280", "Concentration", "A280 → concentration", "")
    R.widget("a280", "a280", "A280 calculator", {"mw": mw, "ext_reduced": red, "ext_cystines": ox, "nW": nW, "nY": nY, "nC": nC},
             how=("Beer–Lambert: c (M) = A280 / (ε × path). ε comes from the Trp, Tyr and cystine counts (Pace et al. 1995; "
                  "Gill & von Hippel 1989). It assumes the protein is pure and that its chromophores absorb as they would unfolded "
                  "(the Pace/Gill–von Hippel approximation, usually within ~5 % for folded proteins). Nucleic acid contamination "
                  "(A260/A280 above ~0.6) inflates A280."),
             yours=(f"Abs 0.1 % = {abs01_red:.3f} (reduced) or {abs01_ox:.3f} (cystines): a 1 mg/mL solution reads A280 ≈ {abs01_red:.2f} in a 1 cm cell."
                    if nW + nY else "No Trp or Tyr: A280 cannot be used for this protein."), wide=False)

    # ---- hydropathy & composition
    R.section("hydro", "Sequence profile", "Hydropathy and composition", "")
    if hyd:
        R.fig("hydropathy", "Kyte–Doolittle hydropathy", "hydropathy", {"windows": {str(k): v for k, v in hyd.items()}, "length": L,
                                                                          "tm": tm, "threshold": 1.6},
              how=("Average Kyte–Doolittle hydropathy in a sliding window (<code>protein_scale(kd, window)</code>). Window 9 shows "
                   "surface-exposed vs buried stretches; window 19 finds transmembrane helices — a 19-residue average above 1.6 "
                   "suggests one (Kyte & Doolittle 1982)."),
              yours=(f"{pl(len(tm), 'candidate transmembrane stretch', 'candidate transmembrane stretches')} (window 19 above 1.6)"
                     + (": " + ", ".join(f"{a}–{b}" for a, b in tm[:6]) if tm else " — consistent with a soluble protein")
                     + f". GRAVY {gravy:+.2f}."), wide=True)
    top = sorted(comp.items(), key=lambda x: -x[1])[:3]
    R.fig("composition", "Amino-acid composition", "bars", {"bars": [[a, round(comp[a], 2)] for a in STD], "unit": "%"},
          how="Percent of each amino acid (<code>amino_acids_percent</code>, with a shim so the browser build's Biopython 1.85 gives the same numbers).",
          yours="Most common: " + ", ".join(f"{a} {v:.1f} %" for a, v in top) + f"; {nC} Cys, {nW} Trp, {seq.count('M')} Met.")

    # ---- tags & protease sites
    rows, found = [], []
    for nm, rx, note in TAGS:
        for m in re.finditer(rx, seq):
            s, e = m.start() + 1, m.end()
            where = "N-terminal" if s <= 25 else "C-terminal" if e >= L - 25 else "internal"
            rows.append([nm, f"{s}–{e}", m.group(), where, note])
            found.append({"kind": "tag", "name": nm, "start": s, "end": e})
    for nm, rx, off, note in PROTEASES:
        for m in re.finditer(rx, seq):
            cut = m.start() + off           # residues before the cut
            down = seq[cut:]
            up = seq[:cut]
            dm = ProteinAnalysis(down).molecular_weight() if len(down) >= 1 else 0
            txt = (f"{note}. After cleavage: upstream {len(up)} aa ending …{up[-6:]}; downstream {len(down)} aa "
                   f"({dm / 1000:.1f} kDa) starting {down[:6]}…")
            rows.append([nm + " site", f"{m.start() + 1}–{m.end()}", m.group(), f"cuts after residue {cut}", txt])
            found.append({"kind": "protease", "name": nm, "start": m.start() + 1, "end": m.end(), "cut": cut})
    R.section("tags", "Tags", "Affinity tags and protease sites", "")
    R.table("tags", "Tags and protease sites", ["What", "Residues", "Sequence", "Where", "Note / scar after cleavage"], rows,
            how=("Exact sequence matches for His6 (six or more H), FLAG, HA, Myc, Strep-tag II and V5, the starts of GST and MBP "
                 "(heuristic: only the first residues are matched), and the TEV, HRV 3C, thrombin and enterokinase cleavage sites. "
                 "The scar is what the protease leaves on the protein downstream of the cut."),
            yours=("None of these tags or protease sites are in this sequence." if not rows else
                   f"{pl(sum(1 for f in found if f['kind'] == 'tag'), 'tag')} and {pl(sum(1 for f in found if f['kind'] == 'protease'), 'protease site')} found."))
    R.fig("seqview", "Sequence", "protseq", {"seq": seq, "marks": found},
          how="Residues numbered from 1; tags are underlined in the accent colour, protease sites in amber with the cut marked ↓.",
          yours=f"{L} residues.", wide=True)
    tags = [f for f in found if f["kind"] == "tag"]
    if tags:
        R.flag("info", "Tag" + ("s" if len(tags) > 1 else "") + " found: " + ", ".join(f"{esc(f['name'])} at {f['start']}–{f['end']}" for f in tags) + ".")
    for f in found:
        if f["kind"] == "protease":
            R.flag("info", f"{esc(f['name'])} site at {f['start']}–{f['end']}: cleavage leaves the downstream protein starting at residue {f['cut'] + 1}.")
    R.method("ProtParam", f"Biopython {__import__('Bio').__version__} <code>ProteinAnalysis</code>: molecular_weight, isoelectric_point "
                          "(Bjellqvist et al. 1993), charge_at_pH, molar_extinction_coefficient (Pace et al. 1995), gravy and "
                          "protein_scale with the Kyte–Doolittle scale (Kyte & Doolittle 1982), instability_index (Guruprasad et al. 1990), "
                          "aromaticity (Lobry & Gautier 1994).")
    R.method("A280", "c = A280 / (ε × l); Abs 0.1 % = ε / MW. Gill & von Hippel 1989; Pace et al. 1995.")
    return R.done()
