"""Primer design: PCR pairs for a region, cloning primers anchored on a CDS, and a sequencing-primer walk.

Built on the same pieces the Primers tool uses to *check* primers — nearest-neighbour Tm under the stated buffer,
the 3′ clamp rule, the complementarity heuristics for dimers and hairpins, and the binding-site search on the whole
template — so a designed primer and a checked primer are judged the same way. Scoring follows primer3's idea
(penalties for distance from the ideal), with the weights written out below; it is not primer3 and computes no
secondary-structure ΔG. Biopython + stdlib only (the browser build runs it too).
"""
from __future__ import annotations

import re

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from .core import Result, UserFacingError, commas, esc, pl
from .core_primers import binding_sites, buffer_from, comp_scan, hairpin_scan, tm_nn
from .core_seq import enzyme, find_feature, record_from_text, record_json, revcomp

MODES = ("pcr", "clone", "seq")
TAIL_PAD = "GCGC"          # extra 5′ bases so the enzyme can cut near the end of the product (NEB: ≥ 4–6 nt)


# ---------------------------------------------------------------- one primer
def assess(seq: str, b: dict, tm_target: float) -> dict:
    """Everything the score needs for one candidate, and why it lost points."""
    s = seq.upper()
    tm = tm_nn(s, b)
    gc = 100 * gc_fraction(s)
    why, pen = [], abs(tm - tm_target)
    if abs(tm - tm_target) >= 2:
        why.append(f"Tm {tm:.1f} °C is {abs(tm - tm_target):.1f} °C from the target")
    if gc < 35 or gc > 65:
        pen += (abs(gc - 50) - 15) * 0.15
        why.append(f"GC {gc:.0f} %")
    n3 = sum(c in "GC" for c in s[-5:])
    if n3 == 0 or n3 >= 4:
        pen += 1.0 if n3 in (0, 4) else 2.0
        why.append("no G/C in the last 5 bases" if n3 == 0 else f"{n3} G/C in the last 5 bases")
    run = max((len(m.group(0)) for m in re.finditer(r"(A+|C+|G+|T+)", s)), default=0)
    if run >= 5:
        pen += 1.5 * (run - 4)
        why.append(f"a run of {run} identical bases")
    di = re.search(r"((..)\2{3,})", s)
    if di:
        pen += 1.5
        why.append(f"dinucleotide repeat {di.group(1)}")
    sd = comp_scan(s, s)
    if sd["three"]["run"] >= 3:
        pen += 1.5 * (sd["three"]["run"] - 2)
        why.append(f"3′ self-complementarity {sd['three']['run']} bp")
    hp = hairpin_scan(s)
    if hp.get("stem", 0) >= 4:
        pen += 1.0 * (hp["stem"] - 3)
        why.append(f"hairpin stem {hp['stem']} bp")
    return {"seq": s, "len": len(s), "tm": round(tm, 2), "gc": round(gc, 1), "pen": round(pen, 2), "why": why,
            "self3": sd["three"]["run"], "hairpin": hp.get("stem", 0)}


def off_targets(seq: str, template: str, circular: bool, intended: tuple) -> list:
    """Other places on the template where the primer's 3′ end would prime (≤ 2 mismatches, last 3 bases exact)."""
    out = []
    for site in binding_sites(seq, template, circular):
        if (site["strand"], site["start0"]) != intended:
            out.append(site)
    return out


# ---------------------------------------------------------------- template and target
def _template(args):
    t = args.get("template") or args.get("record")
    if not t:
        raise UserFacingError("Choose a template: a loaded construct, or paste the sequence.")
    if isinstance(t, str):
        rec, _ = record_from_text(t, bool(args.get("template_circular")), "Pasted template")
        t = record_json(rec)
    return t


def _target(rj: dict, args) -> tuple[int, int, str, dict | None]:
    """(start0, length, label, feature) of the region to work on; start0 is on the original template."""
    L = rj["length"]
    feat = None
    if args.get("feature") not in (None, ""):
        feat = find_feature(rj, args["feature"])
        if not feat:
            raise UserFacingError(f"No feature called “{args['feature']}” in {rj['name']}.")
        s0, e0 = feat["start0"], feat["end0"]
        n = (e0 - s0) % L if feat.get("wraps") else e0 - s0
        return s0, n or L, feat["name"], feat
    try:
        a, b = int(str(args.get("start")).replace(",", "")), int(str(args.get("end")).replace(",", ""))
    except (TypeError, ValueError):
        raise UserFacingError("Give a feature, or a region as start and end (1-based, inclusive).")
    if not (1 <= a <= L and 1 <= b <= L):
        raise UserFacingError(f"The region must lie within 1–{commas(L)}.")
    if b < a and not rj["circular"]:
        raise UserFacingError("On a linear template the end must come after the start.")
    n = (b - a + 1) if b >= a else (L - a + 1 + b)
    return a - 1, n, f"{commas(a)}–{commas(b)}", None


class Frame:
    """The template rotated (if circular) so the target sits at `off` in a linear working string, with enough
    flank on both sides; positions convert back to the original numbering with .orig()."""

    def __init__(self, rj: dict, s0: int, n: int, flank: int):
        self.L, self.circ = rj["length"], rj["circular"]
        seq = rj["seq"].upper()
        if self.circ:
            f = min(flank, max(0, (self.L - n) // 2))
            self.rot = (s0 - f) % self.L
            self.w = (seq[self.rot:] + seq[:self.rot])
            self.off = f
        else:
            self.rot = 0
            self.w = seq
            self.off = s0
        self.n = n

    def orig(self, p: int) -> int:
        return (p + self.rot) % self.L


def _map_feats(rj, frame: Frame, picks: list) -> list:
    """Primer arrows in the map's feature format (original coordinates)."""
    out = []
    for p in picks:
        a = frame.orig(p["w0"])
        b = (a + p["anneal"]) % rj["length"] or rj["length"]
        wraps = b <= a
        out.append({"name": p["name"], "type": "primer_bind", "strand": p["strand"],
                    "parts": [[a, rj["length"]], [0, b]] if wraps else [[a, b]], "start0": a, "end0": b,
                    "wraps": wraps, "length": p["anneal"], "mm": 0, "qualifiers": {}})
    return out


# ---------------------------------------------------------------- PCR
def _candidates(frame: Frame, lo: int, hi: int, strand: int, lens, b, tm_t, keep: int):
    """Primers whose annealing region lies within [lo, hi) of the working string, best first."""
    w, out = frame.w, []
    for ln in lens:
        for i in range(max(0, lo), min(len(w), hi) - ln + 1):
            s = w[i:i + ln] if strand == 1 else revcomp(w[i:i + ln])
            if "N" in s:
                continue
            tm = tm_nn(s, b)
            if abs(tm - tm_t) > 4:
                continue
            out.append((abs(tm - tm_t), i, ln, s))
    out.sort()
    res = []
    for _, i, ln, s in out[:keep * 3]:
        a = assess(s, b, tm_t)
        a.update({"w0": i, "anneal": ln, "strand": strand})
        res.append(a)
    res.sort(key=lambda x: x["pen"])
    return res[:keep]


def design_pcr(R: Result, rj, frame: Frame, label, b, args):
    tm_t = float(args.get("tm") or 60)
    lens = range(int(args.get("min_len") or 18), int(args.get("max_len") or 28) + 1)
    pmin = int(args.get("product_min") or 0) or frame.n + 40
    pmax = int(args.get("product_max") or 0) or frame.n + 400
    if pmax < frame.n + 2 * lens.start:
        raise UserFacingError(f"The largest product ({commas(pmax)} bp) cannot hold the {commas(frame.n)} bp target plus two "
                              "primers; raise the maximum product size.")
    flank = max(40, pmax - frame.n)
    t0, t1 = frame.off, frame.off + frame.n
    fw = _candidates(frame, t0 - flank, t0, 1, lens, b, tm_t, 40)
    rv = _candidates(frame, t1, t1 + flank, -1, lens, b, tm_t, 40)
    if not fw or not rv:
        side = "upstream" if not fw else "downstream"
        raise UserFacingError(f"No primer within 4 °C of {tm_t:.0f} °C fits {side} of the target in the allowed flank. "
                              "Widen the product size range, allow longer primers, or lower the target Tm.")
    pairs = []
    for f in fw:
        for r in rv:
            size = (r["w0"] + r["anneal"]) - f["w0"]
            if not pmin <= size <= pmax:
                continue
            dtm = abs(f["tm"] - r["tm"])
            pairs.append([f["pen"] + r["pen"] + 0.5 * dtm, f, r, size, dtm])
    if not pairs:
        raise UserFacingError(f"No pair gives a product of {commas(pmin)}–{commas(pmax)} bp around the target. Widen the range.")
    pairs.sort(key=lambda x: x[0])
    ranked = []
    for pen, f, r, size, dtm in pairs[:120]:
        cd = comp_scan(f["seq"], r["seq"])
        extra, why = 0.0, []
        if cd["three"]["run"] >= 3:
            extra += 1.5 * (cd["three"]["run"] - 2)
            why.append(f"cross-dimer: {cd['three']['run']} bp at a 3′ end")
        if dtm > 2:
            why.append(f"Tm differ by {dtm:.1f} °C")
        ranked.append([pen + extra, f, r, size, dtm, why])
    ranked.sort(key=lambda x: x[0])
    # off-target priming for the few pairs we show (the costly check, on the whole template both strands)
    out = []
    for pen, f, r, size, dtm, why in ranked:
        if len(out) >= int(args.get("n_pairs") or 5):
            break
        fo = off_targets(f["seq"], rj["seq"].upper(), rj["circular"], (1, frame.orig(f["w0"])))
        ro = off_targets(r["seq"], rj["seq"].upper(), rj["circular"], (-1, frame.orig(r["w0"])))
        if fo or ro:
            pen += 5 * (len(fo) + len(ro))
            why = why + [f"{len(fo) + len(ro)} other binding site(s) on the template"]
        if any(abs(o[1]["w0"] - f["w0"]) < 3 and abs(o[2]["w0"] - r["w0"]) < 3 for o in out):
            continue                                   # near-duplicates of a pair already listed
        if sum(o[1]["seq"] == f["seq"] for o in out) >= 2 or sum(o[2]["seq"] == r["seq"] for o in out) >= 2:
            continue                                   # real alternatives: no primer in more than two listed pairs
        out.append((pen, f, r, size, dtm, why, fo + ro))
    out.sort(key=lambda x: x[0])
    rows, picks = [], []
    for k, (pen, f, r, size, dtm, why, off) in enumerate(out, 1):
        f.update(name=f"F{k}")
        r.update(name=f"R{k}")
        fs, rs = frame.orig(f["w0"]) + 1, frame.orig(r["w0"] + r["anneal"] - 1) + 1
        rows.append([f"Pair {k}", f"{f['seq']}", f"{r['seq']}", f"{commas(fs)} / {commas(rs)}",
                     f"{f['tm']:.1f} / {r['tm']:.1f}", f"{f['gc']:.0f} / {r['gc']:.0f}", commas(size), f"{pen:.1f}",
                     "; ".join(f["why"] + r["why"] + why) or "no issues"])
        if k <= 3:
            picks += [f, r]
    best = out[0]
    R.tile(commas(best[3]) + " bp", "best product")
    R.tile(f"{best[1]['tm']:.1f} / {best[2]['tm']:.1f} °C", "best pair Tm")
    R.tile(len(out), "pairs listed")
    R.r["design"] = {"mode": "pcr", "pairs": [{"fwd": o[1]["seq"], "rev": o[2]["seq"], "size": o[3], "penalty": round(o[0], 2),
                                              "fwd_start0": frame.orig(o[1]["w0"]),
                                              "rev_end0": frame.orig(o[2]["w0"] + o[2]["anneal"] - 1) + 1,
                                              "tm": [o[1]["tm"], o[2]["tm"]]} for o in out],
                     "target": {"start0": frame.orig(frame.off), "length": frame.n}}
    if best[0] < 3:
        R.flag("ok", f"Best pair: product {commas(best[3])} bp, Tm {best[1]['tm']:.1f} and {best[2]['tm']:.1f} °C, "
                     "no complementarity or off-target issues found by the checks below.")
    else:
        R.flag("warn", f"Even the best pair loses points: {esc('; '.join(best[1]['why'] + best[2]['why'] + best[5]) or 'see the table')}. "
                       "Consider a wider product range or a different target Tm.")
    if best[6]:
        R.flag("warn", f"The best pair's primers also match {pl(len(best[6]), 'other site')} on the template (3′ end exact, "
                       "≤ 2 mismatches). The virtual PCR in the Primers tool shows whether that gives an extra product.")
    R.section("pairs", "Primer pairs", f"PCR primers for {esc(label)}",
              f"Candidates within 4 °C of the target Tm on each side of the target, scored, paired for a product of "
              f"{commas(pmin)}–{commas(pmax)} bp, then checked for cross-dimers and for other binding sites on the whole "
              "template (both strands).")
    R.table("pairs", "Best pairs", ["", "Forward 5′→3′", "Reverse 5′→3′", "5′ ends (template)", "Tm F / R (°C)",
                                    "GC F / R (%)", "Product (bp)", "Penalty", "Why it lost points"], rows,
            how="Penalty adds up distances from the ideal: 1 per °C away from the target Tm, 0.5 per °C Tm difference between "
                "the two, and points for GC outside 35–65 %, a 3′ end with 0 or ≥ 4 G/C in the last 5 bases, runs of ≥ 5 "
                "identical bases, dinucleotide repeats, 3′ self-complementarity ≥ 3 bp, hairpin stems ≥ 4 bp, a 3′ cross-dimer "
                "≥ 3 bp, and 5 per other binding site on the template. Lower is better; below about 3 is a clean pair.",
            yours=(f"The best pair amplifies {commas(best[3])} bp including all of {esc(label)} ({commas(frame.n)} bp). "
                   f"Tm {best[1]['tm']:.1f} and {best[2]['tm']:.1f} °C in the buffer used (see Methods) — start Taq annealing near "
                   f"{min(best[1]['tm'], best[2]['tm']) - 4:.0f} °C."))
    return picks, [{"fwd": "F1", "rev": "R1", "start0": frame.orig(best[1]["w0"]),
                    "end0": (frame.orig(best[2]["w0"] + best[2]["anneal"] - 1) + 1) % rj["length"] or rj["length"],
                    "size": best[3], "wraps": False, "mm": 0, "perfect": True}]


# ---------------------------------------------------------------- cloning
def _anchored(frame: Frame, at: int, strand: int, b, tm_t, max_len=40):
    """The shortest primer starting exactly at `at` (forward) or ending exactly at `at` (reverse, exclusive end) that
    reaches the target Tm."""
    best = None
    for ln in range(16, max_len + 1):
        if strand == 1:
            s = frame.w[at:at + ln]
        else:
            if at - ln < 0:
                break
            s = revcomp(frame.w[at - ln:at])
        if len(s) < ln:
            break
        t = tm_nn(s, b)
        best = (s, ln, t)
        if t >= tm_t:
            break
    return best


def design_clone(R: Result, rj, frame: Frame, label, feat, b, args):
    tm_t = float(args.get("tm") or 60)
    t0, t1 = frame.off, frame.off + frame.n
    stop = bool(args.get("keep_stop", True))
    insert = frame.w[t0:t1]
    is_cds = bool(feat and feat["type"] == "CDS")
    notes = []
    if is_cds and feat["strand"] == -1:
        # the CDS is on the bottom strand: the "forward" (start-codon) primer anneals at the right end
        insert = revcomp(insert)
        notes.append("The CDS is on the bottom strand, so the start-codon primer is at the right-hand end of the feature.")
    if is_cds and not stop and insert[-3:] in ("TAA", "TAG", "TGA"):
        insert = insert[:-3]
        notes.append("The stop codon is left out, for a C-terminal fusion.")
    e5, e3 = (args.get("enzyme5") or "").strip(), (args.get("enzyme3") or "").strip()
    tails, checks = {}, []
    kozak = "GCCACC" if args.get("kozak") and is_cds else ""
    for side, name in (("5", e5), ("3", e3)):
        if not name:
            tails[side] = ""
            continue
        e = enzyme(name)
        site = str(e.site)
        if any(c not in "ACGT" for c in site):
            raise UserFacingError(f"{e} recognises a degenerate site ({site}); pick an enzyme with a defined site for a tail.")
        tails[side] = TAIL_PAD + site
        k = site.find("ATG")
        if side == "5" and is_cds and k >= 0 and insert.startswith("ATG"):
            # NcoI (CCATGG), NdeI (CATATG)…: in the vector the site's ATG is the start codon. Putting the whole site in
            # front of the gene's own ATG leaves an ATG out of frame ahead of it, so the two ATGs are merged instead —
            # and any site bases after the ATG overwrite the start of codon 2, which is reported, not hidden.
            rest = site[k + 3:]
            new = insert[:3] + rest + insert[3 + len(rest):]
            if rest and new != insert:
                old_c, new_c = insert[3:6], new[3:6]
                aa_o, aa_n = str(Seq(old_c).translate()), str(Seq(new_c).translate())
                checks.append(("warn", f"<b>{e} fuses its site's ATG with the start codon</b>, and its next base "
                                       f"changes codon 2 from {old_c} ({aa_o}) to {new_c} ({aa_n})"
                                       + (" — a silent change." if aa_o == aa_n else
                                          f" — <b>an amino-acid change {aa_o}2{aa_n}</b>. If residue 2 matters, use "
                                          "NdeI (CATATG), which needs no change after the ATG, or a seamless method.")))
                insert = new
            else:
                checks.append(("ok", f"{e} ({site}) contains the start codon: its ATG is merged with the gene's ATG, so "
                                     "translation starts at the site in frame."))
            tails[side] = TAIL_PAD + site[:k]
            if kozak:
                checks.append(("info", f"No Kozak sequence added: {e}'s site already sets the bases around the ATG."))
                kozak = ""
        cuts = e.search(Seq(insert), linear=True)
        if cuts:
            checks.append(("warn", f"<b>{e} cuts inside the insert</b> ({pl(len(cuts), 'site')}, at "
                                   f"{', '.join(commas(c) for c in cuts[:5])} of the insert): digesting the PCR product "
                                   f"with it would cut the insert too. Pick another enzyme, or use a seamless method "
                                   "(Gibson, In-Fusion)."))
        else:
            checks.append(("ok", f"{e} ({site}) does not cut inside the {commas(len(insert))} bp insert."))
    # the primers are placed after the tails are decided: a fused start-codon site can change the insert's first bases
    ins_frame = Frame({"length": len(insert), "circular": False, "seq": insert}, 0, len(insert), 0)
    f = _anchored(ins_frame, 0, 1, b, tm_t)
    r = _anchored(ins_frame, len(insert), -1, b, tm_t)
    if not f or not r:
        raise UserFacingError("The target is too short to place cloning primers on it.")
    if e5 and e3 and e5.lower() == e3.lower():
        checks.append(("warn", "The same enzyme at both ends: the insert can ligate in either orientation. Screen colonies "
                               "for orientation, or use two different enzymes for directional cloning."))
    fwd = tails["5"].lower() + kozak.lower() + f[0]
    rev = tails["3"].lower() + r[0]
    fa, ra = assess(f[0], b, tm_t), assess(r[0], b, tm_t)
    rows = [["Forward", fwd, f"{f[1]} nt anneal + {len(fwd) - f[1]} nt tail", f"{f[2]:.1f}", f"{fa['gc']:.0f}",
             "; ".join(fa["why"]) or "no issues"],
            ["Reverse", rev, f"{r[1]} nt anneal + {len(rev) - r[1]} nt tail", f"{r[2]:.1f}", f"{ra['gc']:.0f}",
             "; ".join(ra["why"]) or "no issues"]]
    prod = len(insert) + (len(fwd) - f[1]) + (len(rev) - r[1])
    R.tile(commas(prod) + " bp", "PCR product with tails")
    R.tile(f"{f[2]:.1f} / {r[2]:.1f} °C", "Tm of the annealing parts")
    R.r["design"] = {"mode": "clone", "fwd": fwd, "rev": rev, "size": prod}
    for lvl, t in checks:
        R.flag(lvl, t)
    if is_cds and not insert.startswith("ATG"):
        R.flag("warn", f"The insert starts with {insert[:3]}, not ATG — the forward primer does not begin at a start codon.")
    if abs(f[2] - r[2]) > 5:
        R.flag("warn", f"The annealing parts differ by {abs(f[2] - r[2]):.1f} °C in Tm; they are anchored at fixed ends, so "
                       "lengthen the lower one toward the insert or accept a touchdown PCR.")
    for n_ in notes:
        R.flag("info", n_)
    R.section("clone", "Cloning primers", f"Cloning primers for {esc(label)}",
              "Anchored exactly at the two ends of the insert and lengthened base by base until the annealing part reaches the "
              "target Tm. Tails (lower case) are added 5′: extra bases so the enzyme can cut near the end, then the site"
              + (", then a Kozak sequence (GCCACC) before the ATG" if kozak else "") + ".")
    R.table("cloneprimers", "Primers", ["", "Sequence 5′→3′ (tail lower case)", "Parts", "Tm anneal (°C)", "GC (%)", "Notes"],
            rows, how="The Tm is for the part that anneals in the first cycles; once the tails are copied into the product, "
                      "later cycles anneal over the full length, so a 2-step protocol (a few cycles at the lower Tm, then "
                      "higher) is common.",
            yours=f"Product {commas(prod)} bp: the {commas(len(insert))} bp insert plus {len(fwd) - f[1]} + {len(rev) - r[1]} nt of tails.")
    fr = {"name": "F", "seq": fwd, "w0": frame.off if not (is_cds and feat["strand"] == -1) else t1 - f[1],
          "anneal": f[1], "strand": 1 if not (is_cds and feat["strand"] == -1) else -1}
    rr = {"name": "R", "seq": rev, "w0": t1 - r[1] if not (is_cds and feat["strand"] == -1) else t0,
          "anneal": r[1], "strand": -1 if not (is_cds and feat["strand"] == -1) else 1}
    return [fr, rr], []


# ---------------------------------------------------------------- sequencing walk
def design_seq(R: Result, rj, frame: Frame, label, b, args):
    read = int(args.get("read_len") or 700)
    if read < 200:
        raise UserFacingError("A usable read length below 200 bp is not a Sanger read; use 600–800.")
    tm_t = float(args.get("tm") or 56)
    lens = range(18, 25)
    both = bool(args.get("both_strands"))
    dead = 50                                    # the first ~50 bases after a primer are unreadable
    step = read - 100                            # 100 bp overlap between neighbouring reads
    picks, reads = [], []
    t0, t1 = frame.off, frame.off + frame.n
    for strand in ((1, -1) if both else (1,)):
        k, cover = 0, (t0 if strand == 1 else t1)
        while (cover < t1) if strand == 1 else (cover > t0):
            k += 1
            if strand == 1:
                want = cover - dead - 20           # primer 3′ end ~dead bp before what this read must cover
                cands = _candidates(frame, want - 60, want + 20, 1, lens, b, tm_t, 5)
            else:
                want = cover + dead
                cands = _candidates(frame, want - 20, want + 60, -1, lens, b, tm_t, 5)
            if not cands:
                raise UserFacingError(f"No sequencing primer near {tm_t:.0f} °C fits around position "
                                      f"{commas(frame.orig(max(0, want)) + 1)}; lower the target Tm.")
            c = cands[0]
            c["name"] = f"{'S' if strand == 1 else 'A'}{k}"
            picks.append(c)
            if strand == 1:
                s_ = c["w0"] + c["anneal"] + dead
                e_ = min(len(frame.w), c["w0"] + c["anneal"] + read)
                reads.append((c, s_, e_))
                cover = max(cover + 50, e_ - 100)
            else:
                e_ = c["w0"] - dead
                s_ = max(0, c["w0"] - read)
                reads.append((c, s_, e_))
                cover = min(cover - 50, s_ + 100)
            if k > 60:
                break
    rows = []
    for c, s_, e_ in reads:
        rows.append([c["name"], "forward" if c["strand"] == 1 else "reverse", c["seq"],
                     commas(frame.orig(c["w0"]) + 1), f"{c['tm']:.1f}", f"{c['gc']:.0f}",
                     f"{commas(frame.orig(s_) + 1)}–{commas(frame.orig(e_ - 1) + 1)}", "; ".join(c["why"]) or "no issues"])
    n_f = sum(1 for c in picks if c["strand"] == 1)
    R.tile(len(picks), "sequencing primers")
    R.tile(commas(frame.n) + " bp", "region to cover")
    R.tile("both strands" if both else "one strand", "coverage")
    R.r["design"] = {"mode": "seq", "primers": [{"name": c["name"], "seq": c["seq"], "strand": c["strand"],
                                                 "start0": frame.orig(c["w0"])} for c in picks],
                     "reads": [{"primer": c["name"], "strand": c["strand"], "from0": frame.orig(s_), "len": e_ - s_}
                               for c, s_, e_ in reads],
                     "target": {"start0": frame.orig(frame.off), "length": frame.n}}
    R.flag("ok", f"{pl(n_f, 'forward primer')}" + (f" and {pl(len(picks) - n_f, 'reverse primer')}" if both else "")
                 + f" cover {esc(label)} ({commas(frame.n)} bp) with ~{read} bp reads overlapping by ~100 bp.")
    if not both:
        R.flag("info", "One strand only: a difference seen on one read is not confirmed by an independent read. For a "
                       "construct you will publish or ship, sequence both strands (tick “both strands”).")
    R.section("seq", "Sequencing", f"Sequencing primers across {esc(label)}",
              f"A primer every ~{step} bp: each sits ~{dead + 20} bp before the stretch its read must cover (the first ~{dead} "
              f"bases of a Sanger read are unreadable) and reads ~{read} bp.")
    R.table("seqprimers", "Primers and the stretch each read covers",
            ["", "Strand", "Sequence 5′→3′", "5′ end", "Tm (°C)", "GC (%)", "Expected good read", "Notes"], rows,
            how="Expected reads assume ~" + str(read) + " good bases per reaction; real reads vary with template quality.",
            yours=f"{pl(len(picks), 'primer')} for {commas(frame.n)} bp.")
    prods = [{"fwd": c["name"], "rev": c["name"], "start0": frame.orig(s_), "end0": frame.orig(e_ - 1) + 1, "size": e_ - s_,
              "wraps": frame.orig(e_ - 1) < frame.orig(s_), "mm": 0, "perfect": True} for c, s_, e_ in reads]
    return picks, prods


# ---------------------------------------------------------------- the tool
def design(args: dict) -> dict:
    mode = (args.get("mode") or "pcr").lower()
    if mode not in MODES:
        raise UserFacingError(f"Unknown design mode “{mode}” (use pcr, clone or seq).")
    rj = _template(args)
    s0, n, label, feat = _target(rj, args)
    if mode == "clone" and n > 20000:
        raise UserFacingError("That insert is over 20 kb — too long for one PCR.")
    b = buffer_from(args)
    flank = int(args.get("product_max") or 0) or n + 400
    frame = Frame(rj, s0, n, max(400, flank))
    R = Result("design", f"{ {'pcr': 'PCR primers', 'clone': 'Cloning primers', 'seq': 'Sequencing primers'}[mode]} · {label}")
    R.r["record"] = rj
    R.r["buffer"] = b
    if mode == "pcr":
        picks, prods = design_pcr(R, rj, frame, label, b, args)
    elif mode == "clone":
        picks, prods = design_clone(R, rj, frame, label, feat, b, args)
    else:
        picks, prods = design_seq(R, rj, frame, label, b, args)
    R.section("map", "Map", f"Where the primers sit on {esc(rj['name'])}", "")
    R.fig("designmap", "Primers on the template", "map", {"record_ref": True, "primers": _map_feats(rj, frame, picks),
                                                          "products": prods, "cutters": []},
          how="Arrows are the annealing parts of the primers, pointing in the direction of extension"
              + ("; the shaded arcs are the reads each sequencing primer is expected to give." if mode == "seq"
                 else "; the shaded arc is the product of the best pair." if prods else "."),
          yours=f"Target {esc(label)} ({commas(n)} bp) on the {'circular' if rj['circular'] else 'linear'} "
                f"{commas(rj['length'])} bp template.")
    R.method("Primer design", "Candidates scored with nearest-neighbour Tm (Biopython <code>MeltingTemp.Tm_NN</code>, DNA_NN3 "
                              "Allawi &amp; SantaLucia 1997; Owczarzy 2008 salt correction with Mg²⁺ and dNTPs) under the buffer "
                              f"{esc(b['label'])} (K⁺ {b['K']:g}, Na⁺ {b['Na']:g}, Mg²⁺ {b['Mg']:g}, dNTPs {b['dNTPs']:g} mM; primer "
                              f"{b['primer_nM']:g} nM), GC, 3′ clamp, repeats, and a complementarity heuristic for dimers and "
                              "hairpins (no ΔG); penalties in the style of primer3 (Untergasser et al. 2012) with the weights "
                              "stated in the table. Off-target sites: ≤ 2 mismatches with the last 3 3′ bases exact, both strands.")
    R.r["summary"] = {"mode": mode, "target": label, "template": rj["name"]}
    return R.done()
