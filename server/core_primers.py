"""Primers: Tm by three methods under a stated buffer, 3′ clamp, repeats, complementarity heuristics, pair
analysis and a virtual PCR against a template. Biopython + stdlib only.

Biopython has no nucleic-acid secondary-structure thermodynamics (no ΔG for dimers or hairpins). The
dimer and hairpin checks here are a transparent complementarity heuristic and say so everywhere they appear.
"""
from __future__ import annotations

import itertools
import re

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction, molecular_weight

from .core import Result, UserFacingError, commas, esc, pl
from .core_seq import features_at, normalize_dna, revcomp

PRESETS = {
    "taq": {"label": "Taq, standard", "Na": 0, "K": 50, "Tris": 10, "Mg": 1.5, "dNTPs": 0.8, "primer_nM": 250, "DMSO": 0,
            "why": ("A standard Taq buffer: 50 mM KCl, 10 mM Tris-HCl, 1.5 mM MgCl₂, 0.2 mM of each dNTP (0.8 mM total) and "
                    "250 nM of each primer. dNTPs chelate Mg²⁺, so the free Mg²⁺ that stabilises the duplex is lower than 1.5 mM; "
                    "the Owczarzy 2008 correction accounts for that.")},
    "hifi": {"label": "High-fidelity (Q5/Phusion-like)", "Na": 0, "K": 0, "Tris": 0, "Mg": 2.0, "dNTPs": 0.8, "primer_nM": 500, "DMSO": 0,
             "why": ("2 mM Mg²⁺, 0.2 mM of each dNTP, 500 nM of each primer, as in Q5 and Phusion reactions. NEB's and Thermo's own "
                     "Tm calculators use their own proprietary parameter sets and buffer compositions, so the values here will "
                     "differ from theirs by a few °C. For Q5 or Phusion annealing temperatures, use the vendor calculator "
                     "(tmcalculator.neb.com, or Thermo's Tm calculator) — they also add the polymerase's own rule (e.g. Tm + 3 °C).")},
}
IUPAC_EXPAND = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
                "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
PAIR = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}


# ---------------------------------------------------------------- parsing
def parse_primers(text: str) -> list:
    out = []
    for n, line in enumerate((text or "").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            out.append({"name": line[1:].strip() or f"P{len(out) + 1}", "raw": ""})
            continue
        parts = re.split(r"[\t,;]+|\s{1,}", line)
        # A first word is a name unless it could only be sequence. "R ACGT…" (R is also an IUPAC code) and "KanR TTAG…"
        # (every letter of KanR is one) are how people name primers; read as sequence they glued extra "bases" onto
        # the 5′ end. Sequence is written in one case, so a mixed-case word, or one of three letters or fewer ahead of
        # a real sequence, is a name.
        first = parts[0]
        named = len(parts) >= 2 and len("".join(parts[1:])) >= 10 and (
            len(first) <= 3 or not (first.isupper() or first.islower()))
        if named or (len(parts) >= 2 and not re.fullmatch(r"[ACGTUacgtuRYSWKMBDHVNryswkmbdhvn\-]+", first)):
            name, seq = parts[0], "".join(parts[1:])
        elif len(parts) >= 2 and all(re.fullmatch(r"[ACGTUacgtuRYSWKMBDHVNryswkmbdhvn]+", p) for p in parts):
            name, seq = None, "".join(parts)
        else:
            name, seq = (parts[0], "".join(parts[1:])) if len(parts) > 1 else (None, parts[0])
        seq = re.sub(r"^5'?-?|-?3'?$", "", seq.replace("5′", "").replace("3′", ""))
        if out and out[-1]["raw"] == "" and name is None:
            out[-1]["raw"] = seq
            continue
        out.append({"name": name or f"P{len(out) + 1}", "raw": seq})
    out = [p for p in out if p["raw"]]
    if not out:
        raise UserFacingError("Paste at least one primer (one per line, optionally “name sequence”).")
    if len(out) > 40:
        raise UserFacingError("Paste at most 40 primers at a time.")
    names = set()
    for p in out:
        base, k = p["name"], 2
        while p["name"] in names:
            p["name"] = f"{base}_{k}"
            k += 1
        names.add(p["name"])
    return out


def buffer_from(args: dict) -> dict:
    """The reaction buffer: a preset, with any values the user changed laid over it (then labelled Custom)."""
    preset = (args.get("preset") or "taq").lower()
    base = PRESETS.get(preset, PRESETS["taq"])
    b = dict(base)
    user = args.get("buffer") or {}
    changed = False
    for k in ("Na", "K", "Tris", "Mg", "dNTPs", "primer_nM", "DMSO"):
        if user.get(k) in (None, ""):
            continue
        try:
            v = float(user[k])
        except (TypeError, ValueError):
            raise UserFacingError(f"The buffer value for {k} is not a number.")
        if v < 0:
            raise UserFacingError(f"The buffer value for {k} cannot be negative.")
        if abs(v - float(base[k])) > 1e-9:
            changed = True
        b[k] = v
    if preset == "custom" or changed:
        b["label"] = "Custom" if preset not in PRESETS else f"Custom (from {base['label']})"
        b["why"] = "Your own buffer values." + (" " + base["why"] if preset == "hifi" else "")
    if b["Na"] + b["K"] + b["Tris"] + b["Mg"] <= 0:
        raise UserFacingError("The buffer has no salt at all (Na⁺, K⁺, Tris and Mg²⁺ are all 0). Tm is undefined without cations.")
    if b["Na"] + b["K"] + b["Tris"] <= 0 and b["Mg"] <= b["dNTPs"]:
        # dNTPs bind Mg²⁺ one to one, so with no monovalent salt nothing is left to stabilise the duplex and the
        # salt correction divides by zero (Biopython raised "Total ion concentration of zero")
        raise UserFacingError(f"All the Mg²⁺ ({b['Mg']:g} mM) is taken up by the dNTPs ({b['dNTPs']:g} mM total), and there is "
                              "no Na⁺, K⁺ or Tris — so no free cations and no defined Tm. PCR needs Mg²⁺ above the total "
                              "dNTP concentration (typically 1.5–2 mM Mg²⁺ with 0.8 mM total dNTPs).")
    if b["primer_nM"] <= 0:
        raise UserFacingError("Primer concentration must be above 0 nM.")
    if b["DMSO"] > 20:
        raise UserFacingError("DMSO above 20 % is outside any PCR buffer; check the value.")
    return b


# ---------------------------------------------------------------- per-primer numbers
def tm_nn(seq: str, b: dict) -> float:
    """Nearest-neighbour Tm: Allawi & SantaLucia 1997 DNA_NN3, Owczarzy 2008 salt correction (saltcorr=7, handles
    Mg²⁺ and dNTP chelation), primer in excess over template (dnac1 = primer nM, dnac2 = 0), then DMSO."""
    t = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, saltcorr=7, Na=b["Na"], K=b["K"], Tris=b["Tris"], Mg=b["Mg"],
                 dNTPs=b["dNTPs"], dnac1=b["primer_nM"], dnac2=0)
    if b.get("DMSO"):
        t = mt.chem_correction(t, DMSO=b["DMSO"])
    return t


def tm_gc(seq: str, b: dict) -> float:
    t = mt.Tm_GC(seq, Na=b["Na"], K=b["K"], Tris=b["Tris"], Mg=b["Mg"], dNTPs=b["dNTPs"], valueset=7)
    if b.get("DMSO"):
        t = mt.chem_correction(t, DMSO=b["DMSO"])
    return t


def expand(seq: str, cap: int = 256) -> list:
    pools = [IUPAC_EXPAND[c] for c in seq]
    n = 1
    for p in pools:
        n *= len(p)
    if n > cap:
        return []
    return ["".join(x) for x in itertools.product(*pools)]


def runs_and_repeats(s: str) -> list:
    out = []
    for m in re.finditer(r"(A{4,}|C{4,}|G{4,}|T{4,})", s):
        out.append(f"run of {len(m.group())} {m.group()[0]} at {m.start() + 1}")
    for m in re.finditer(r"((?:[ACGT]{2}))\1{3,}", s):
        unit = m.group(1)
        if unit[0] != unit[1]:
            out.append(f"{len(m.group()) // 2}× {unit} repeat at {m.start() + 1}")
    return out


def clamp(s: str) -> tuple[str, str]:
    last5 = s[-5:]
    n = sum(1 for c in last5 if c in "GC")
    if n == 0:
        return "info", f"no G/C in the last 5 bases ({last5}): a weak 3′ end primes less efficiently"
    if n <= 2:
        return "ok", f"{n} G/C in the last 5 bases ({last5}): ideal"
    if n == 3:
        return "ok", f"3 G/C in the last 5 bases ({last5}): acceptable"
    return "warn", f"{n} G/C in the last 5 bases ({last5}): a very GC-rich 3′ end can prime at partial matches"


# ---------------------------------------------------------------- complementarity heuristics
def comp_scan(a: str, b: str) -> dict:
    """Slide primer b (antiparallel) along primer a and find (1) the longest run of Watson–Crick pairs and
    (2) the longest run that starts at either 3′ end. A heuristic: counts pairs, computes no free energy."""
    br = b[::-1]            # b written 3'→5', so pairing positions line up with a written 5'→3'
    la, lb = len(a), len(br)
    best = {"run": 0}
    best3 = {"run": 0}
    for shift in range(-(lb - 1), la):
        run = 0
        for i in range(max(0, shift), min(la, shift + lb)):
            j = i - shift
            if (a[i], br[j]) in PAIR:
                run += 1
                if run > best["run"]:
                    best = {"run": run, "shift": shift, "end_i": i}
            else:
                run = 0
        # 3′-anchored runs: a's 3′ end (i = la-1) or b's 3′ end (j = 0, i = shift)
        for anchor in ("a", "b"):
            if anchor == "a":
                i, j = la - 1, la - 1 - shift
                if not 0 <= j < lb:
                    continue
                k = 0
                while i - k >= 0 and j - k >= 0 and (a[i - k], br[j - k]) in PAIR:
                    k += 1
                if k > best3["run"]:
                    best3 = {"run": k, "shift": shift, "end_i": la - 1, "anchor": "a"}
            else:
                i, j = shift, 0
                if not 0 <= i < la:
                    continue
                k = 0
                while i + k < la and j + k < lb and (a[i + k], br[j + k]) in PAIR:
                    k += 1
                if k > best3["run"]:
                    best3 = {"run": k, "shift": shift, "end_i": shift + k - 1, "anchor": "b"}
    return {"longest": best, "three": best3, "art": dimer_art(a, b, (best3 if best3["run"] >= 3 else best)) if best["run"] else ""}


def dimer_art(a: str, b: str, hit: dict) -> str:
    """Two-line picture of b pairing antiparallel under a, bars at paired bases of the reported run."""
    br = b[::-1]
    shift = hit["shift"]
    pad_a = max(0, -shift)
    pad_b = max(0, shift)
    top = " " * pad_a + a
    bot = " " * pad_b + br
    run = hit["run"]
    end_i = hit["end_i"] + pad_a
    mid = [" "] * max(len(top), len(bot))
    for x in range(end_i - run + 1, end_i + 1):
        if 0 <= x < len(mid):
            mid[x] = "|"
    for x in range(min(len(top), len(bot))):
        if mid[x] == " " and top[x] != " " and bot[x] != " " and (top[x], bot[x]) in PAIR:
            mid[x] = "·"
    return f"5′ {top} 3′\n   {''.join(mid).rstrip()}\n3′ {bot} 5′"


def hairpin_scan(s: str, min_loop: int = 3) -> dict:
    """Longest stem (consecutive intramolecular pairs) with a loop of at least min_loop bases."""
    n = len(s)
    best = {"stem": 0}
    for i in range(n):
        for j in range(n - 1, i + min_loop, -1):
            k = 0
            while i + k < j - k - min_loop and (s[i + k], s[j - k]) in PAIR:
                k += 1
            if k > best["stem"]:
                loop = (j - k) - (i + k) + 1
                best = {"stem": k, "i": i, "j": j, "loop": loop, "three_prime": j == n - 1}
    if best["stem"]:
        i, j, k = best["i"], best["j"], best["stem"]
        arm5 = s[i:i + k]
        arm3 = s[j - k + 1:j + 1][::-1]
        loop = s[i + k:j - k + 1]
        pre = s[:i]
        post = s[j + 1:]
        line1 = f"5′ {pre}{arm5}"
        line3 = f"3′ {post[::-1]}{arm3}"
        w = max(len(line1), len(line3))
        line1 = line1.rjust(w) + "─┐"
        line3 = line3.rjust(w) + "─┘"
        mid = " " * (w - k) + "|" * k + f"  loop {loop}"
        best["art"] = f"{line1}\n{mid}\n{line3}"
    return best


def primer_row(p: dict, b: dict) -> dict:
    raw = p["raw"]
    mixed = any(c.islower() for c in raw) and any(c.isupper() for c in raw)
    tail = 0
    if mixed:
        m = re.match(r"^([a-z]+)([A-Z].*)$", raw)
        tail = len(m.group(1)) if m else 0
    s, notes = normalize_dna(raw, f"Primer {p['name']}")
    r = {"name": p["name"], "seq": s, "length": len(s), "flags": [], "tail": tail}
    for lvl, t in notes:
        if "lower-case" in t:
            continue
        r["flags"].append((lvl, re.sub(r"<[^>]+>", "", t)))
    amb = [c for c in s if c not in "ACGT"]
    r["gc"] = 100 * gc_fraction(s, ambiguous="ignore")
    if amb:
        r["flags"].append(("warn", f"contains IUPAC code{'s' if len(amb) > 1 else ''} {''.join(sorted(set(amb)))} — a degenerate primer is a "
                                   "mix; Tm is given as the range over its variants"))
        vs = expand(s)
        if vs:
            tms = [tm_nn(v, b) for v in vs]
            r["tm_nn"], r["tm_nn_range"] = sum(tms) / len(tms), (min(tms), max(tms))
            gcs = [tm_gc(v, b) for v in vs]
            r["tm_gc"] = sum(gcs) / len(gcs)
            r["tm_w"] = sum(mt.Tm_Wallace(v) for v in vs) / len(vs)
            r["mw"] = sum(molecular_weight(Seq(v), "DNA") for v in vs) / len(vs)
        else:
            r["tm_nn"] = r["tm_gc"] = r["tm_w"] = r["mw"] = None
            r["flags"].append(("warn", "too degenerate to enumerate (> 256 variants): Tm not computed"))
    else:
        r["tm_nn"] = tm_nn(s, b)
        r["tm_gc"] = tm_gc(s, b)
        r["tm_w"] = mt.Tm_Wallace(s)
        r["mw"] = molecular_weight(Seq(s), "DNA")
    if tail:
        bind = s[tail:]
        r["tm_bind"] = tm_nn(bind, b) if not set(bind) - set("ACGT") and len(bind) >= 8 else None
        r["flags"].append(("info", f"lower-case 5′ part ({tail} nt) read as a tail; the upper-case {len(bind)} nt that anneal in "
                                   f"the first cycles have Tm {r['tm_bind']:.1f} °C" if r.get("tm_bind") else
                                   f"lower-case 5′ part ({tail} nt) read as a tail"))
    if len(s) < 15:
        r["flags"].append(("warn", f"only {len(s)} nt: short primers bind at many places"))
    elif len(s) > 60:
        r["flags"].append(("info", f"{len(s)} nt: long oligos are more likely to carry synthesis errors; check the purification"))
    if r["gc"] < 35 or r["gc"] > 65:
        r["flags"].append(("warn", f"GC {r['gc']:.0f} % is outside the usual 40–60 %"))
    lvl, t = clamp(s)
    r["clamp"] = t
    r["clamp_level"] = lvl
    if lvl != "ok":
        r["flags"].append((lvl, t))
    rr = runs_and_repeats(s)
    r["repeats"] = rr
    for x in rr:
        severe = not x.startswith("run of 4 ")
        r["flags"].append(("warn" if severe else "info", x + (" (slippage / mispriming risk)" if severe else " (a mild slippage risk)")))
    # complementarity
    sd = comp_scan(s, s)
    r["self"] = sd
    if sd["three"]["run"] >= 3:
        r["flags"].append(("warn", f"self-dimer: its 3′ end pairs with itself over {sd['three']['run']} bp — can extend into primer-dimer"))
    elif sd["longest"]["run"] >= 8:
        r["flags"].append(("warn", f"self-complementary stretch of {sd['longest']['run']} bp"))
    hp = hairpin_scan(s)
    r["hairpin"] = hp
    if hp["stem"] >= 4:
        r["flags"].append(("warn" if hp["three_prime"] or hp["stem"] >= 6 else "info",
                           f"possible hairpin: {hp['stem']} bp stem, {hp['loop']} nt loop"
                           + (", 3′ end in the stem (it can self-prime)" if hp["three_prime"] else "")))
    return r


# ---------------------------------------------------------------- PCR against a template
def binding_sites(primer: str, template: str, circular: bool, max_mm: int = 2, exact3: int = 3) -> list:
    """Where a primer anneals: ≤ max_mm mismatches over its length, and its last exact3 3′ bases matching
    exactly. Forward sites anneal to the bottom strand and extend rightwards; reverse sites extend leftwards."""
    n, L = len(primer), len(template)
    if n > L or n < 8:
        return []
    ext = template + template[:n - 1] if circular else template
    rc = revcomp(primer)
    out = []
    for strand, probe in ((1, primer), (-1, rc)):
        for pos in range(0, (L if circular else L - n + 1)):
            win = ext[pos:pos + n]
            if len(win) < n:
                break
            # 3′ end first: cheap rejection
            if strand == 1:
                if win[-exact3:] != probe[-exact3:]:
                    continue
            else:
                if win[:exact3] != probe[:exact3]:
                    continue
            mm = 0
            for x, y in zip(win, probe):
                if x != y and not (y in IUPAC_EXPAND and x in IUPAC_EXPAND[y]):
                    mm += 1
                    if mm > max_mm:
                        break
            if mm <= max_mm:
                out.append({"strand": strand, "start0": pos, "end0": (pos + n - 1) % L + 1, "mm": mm})
    return out


def pcr(primers: list, rj: dict, max_product: int = 15000) -> dict:
    seq, L, circ = rj["seq"], rj["length"], rj["circular"]
    sites = []
    for p in primers:
        s = p["seq"]
        found = binding_sites(s, seq, circ)
        tail = 0
        if not found and len(s) > 22:
            # a cloning primer: only its 3′ part matches the template; the 5′ part is a tail. Place it by its 3′ 20 nt,
            # then extend the match 5′-ward base by base to find where the tail really begins.
            core3 = s[-20:]
            found = binding_sites(core3, seq, circ)
            ext = seq + seq if circ else seq
            for f in found:
                k = 20
                if f["strand"] == 1:
                    while k < len(s) and (f["start0"] - 1 >= 0 or circ) and ext[(f["start0"] - 1) % L] == s[-k - 1]:
                        f["start0"] = (f["start0"] - 1) % L
                        k += 1
                else:
                    rcp = revcomp(s)           # its 3′ end is rcp's start; the tail extends past the site's right end
                    while k < len(s) and (f["end0"] < L or circ) and ext[f["end0"] % L] == rcp[k]:
                        f["end0"] = f["end0"] % L + 1
                        k += 1
                f["tail_len"] = len(s) - k
            tail = None
        elif p.get("tail"):
            core3 = s[p["tail"]:]
            f2 = binding_sites(core3, seq, circ)
            if f2 and not found:
                found, tail = f2, p["tail"]
        for f in found:
            f.update({"primer": p["name"], "tail": f.pop("tail_len", tail or 0), "len": len(s)})
            f["features"] = features_at(rj, f["start0"] + 1) or features_at(rj, (f["end0"] - 1) % L)
        sites += found
    fwd = [x for x in sites if x["strand"] == 1]
    rev = [x for x in sites if x["strand"] == -1]
    products = []
    for f in fwd:
        for r in rev:
            a = f["start0"]
            b = r["end0"]                      # exclusive end on the top strand
            if circ:
                size = (b - a) % L or L
                if b <= a:
                    size = L - a + b
            else:
                if b <= a:
                    continue
                size = b - a
            size_tail = size + f["tail"] + r["tail"]
            if size < f["len"] - f["tail"] or size > max_product:
                continue
            products.append({"fwd": f["primer"], "rev": r["primer"], "start0": a, "end0": b, "size": size_tail,
                             "template_span": size, "wraps": circ and b <= a, "mm": f["mm"] + r["mm"],
                             "perfect": f["mm"] == 0 and r["mm"] == 0})
    products.sort(key=lambda p: (p["mm"], p["size"]))
    return {"sites": sites, "products": products}


# ---------------------------------------------------------------- the tool
def primers(args: dict) -> dict:
    plist = parse_primers(args.get("text") or "")
    b = buffer_from(args)
    rows = [primer_row(p, b) for p in plist]
    name = ", ".join(r["name"] for r in rows[:3]) + (f" +{len(rows) - 3}" if len(rows) > 3 else "")
    R = Result("primers", name)
    R.r["primers"] = [{k: v for k, v in r.items() if k in ("name", "seq", "length", "gc", "tm_nn", "tm_gc", "tm_w", "mw", "tail")} for r in rows]
    R.r["buffer"] = b
    tms = [r["tm_nn"] for r in rows if r.get("tm_nn") is not None]
    R.tile(len(rows), "primers")
    if tms:
        R.tile(f"{min(tms):.1f}–{max(tms):.1f} °C" if len(tms) > 1 else f"{tms[0]:.1f} °C", "Tm (nearest-neighbour)")
    R.tile(b["label"], "buffer")
    n_warn = sum(1 for r in rows if any(l == "warn" for l, _ in r["flags"]))
    R.tile(n_warn, "primers with a CHECK")
    if b["Mg"] <= b["dNTPs"]:
        R.flag("warn", f"The buffer has {b['Mg']:g} mM Mg²⁺ and {b['dNTPs']:g} mM total dNTPs. dNTPs bind Mg²⁺ one to one, so "
                       "there is no free Mg²⁺ left for the polymerase and the PCR is unlikely to work; the Tm shown assumes "
                       "the monovalent salt alone. Is the dNTP value per nucleotide rather than total (0.2 mM each = 0.8 mM)?")
    for r in rows:
        for lvl, t in r["flags"]:
            if lvl in ("warn", "info"):
                R.flag(lvl, f"<b>{esc(r['name'])}</b>: {esc(t)}.")
    clean = [r["name"] for r in rows if not r["flags"]]
    if clean:
        R.flag("ok", f"{', '.join(esc(c) for c in clean)}: no issues found by the length, GC, clamp, repeat and complementarity checks.")

    R.section("props", "Per primer", "Melting temperature and composition",
              f"Buffer: {esc(b['label'])} — Na⁺ {b['Na']:g} mM, K⁺ {b['K']:g} mM, Tris {b['Tris']:g} mM, Mg²⁺ {b['Mg']:g} mM, "
              f"dNTPs {b['dNTPs']:g} mM, primer {b['primer_nM']:g} nM, DMSO {b['DMSO']:g} %.")
    trows = []
    for r in rows:
        tmn = "—" if r["tm_nn"] is None else (f"{r['tm_nn']:.1f}" + (f" ({r['tm_nn_range'][0]:.1f}–{r['tm_nn_range'][1]:.1f})" if r.get("tm_nn_range") else ""))
        trows.append([r["name"], r["seq"], r["length"], f"{r['gc']:.1f}", tmn,
                      "—" if r["tm_gc"] is None else f"{r['tm_gc']:.1f}", "—" if r["tm_w"] is None else f"{r['tm_w']:.0f}",
                      "—" if r["mw"] is None else f"{r['mw']:,.1f}", r["clamp"].split(":")[0],
                      "; ".join(r["repeats"]) or "none"])
    diff = (max(tms) - min(tms)) if len(tms) > 1 else 0
    R.table("primers", "Primers", ["Name", "Sequence 5′→3′", "nt", "GC %", "Tm NN (°C)", "Tm GC (°C)", "Tm Wallace (°C)",
                                    "MW (g/mol, ssDNA)", "3′ clamp", "Runs / repeats"], trows,
            how=("<b>Tm NN</b>: nearest-neighbour thermodynamics (<code>Tm_NN</code>, table DNA_NN3 = Allawi & SantaLucia 1997) "
                 "with the Owczarzy 2008 salt correction (<code>saltcorr=7</code>, which models Mg²⁺ and its chelation by dNTPs), "
                 "primer in excess over template (dnac1 = primer nM, dnac2 = 0), then <code>chem_correction</code> for DMSO "
                 "(−0.75 °C per %). <b>Tm GC</b>: an empirical GC-content formula (<code>Tm_GC</code>, valueset 7, the Primer3Plus "
                 "variant, with a sodium-equivalent salt term). <b>Tm Wallace</b>: 2 °C per A/T + 4 °C per G/C. "
                 "<b>Why they differ</b>: only the NN method accounts for which bases neighbour which, the salts and the primer "
                 "concentration; the GC formula knows only composition and length, and the Wallace rule only counts bases — "
                 "it is meant for 14–20-mers in 0.9 M NaCl and overestimates longer primers."),
            yours=(f"Nearest-neighbour Tm ranges {min(tms):.1f}–{max(tms):.1f} °C across {pl(len(tms), 'primer')}"
                   + (f" (a spread of {diff:.1f} °C)" if len(tms) > 1 else "") + "." if tms else "No Tm could be computed."))

    R.section("structure", "Complementarity", "Self-dimers, hairpins and cross-dimers",
              "<b>A heuristic, not thermodynamics.</b> Biopython has no secondary-structure ΔG. This check slides each primer "
              "against itself and its partners and reports the longest run of Watson–Crick pairs, with extra weight on runs "
              "that include a 3′ end (those can be extended by the polymerase); hairpins are stems of ≥ 4 bp closing a loop of "
              "≥ 3 nt. primer3 (ntthal) or the vendor tools compute ΔG — use them to decide whether a flagged structure is stable "
              "at your annealing temperature.")
    for r in rows:
        sd, hp = r["self"], r["hairpin"]
        txt = []
        if sd["longest"]["run"]:
            txt.append(f"Self-dimer (longest run {sd['longest']['run']} bp; 3′-anchored run {sd['three']['run']} bp):\n{sd['art']}")
        if hp.get("stem", 0) >= 3:
            txt.append(f"Hairpin (stem {hp['stem']} bp, loop {hp['loop']} nt{', 3′ end in the stem' if hp['three_prime'] else ''}):\n{hp['art']}")
        verdict = []
        verdict.append(f"3′ self-complementarity {sd['three']['run']} bp — " + ("<b>CHECK</b>: ≥ 3 bp at the 3′ end can form primer-dimer."
                                                                             if sd["three"]["run"] >= 3 else "fine (< 3 bp)."))
        verdict.append(f"Longest self-complementary run {sd['longest']['run']} bp.")
        verdict.append(f"Longest hairpin stem {hp.get('stem', 0)} bp" + (" (≥ 4 bp with a ≥ 3 nt loop: worth a ΔG check)." if hp.get("stem", 0) >= 4 else " (below the 4 bp threshold)."))
        R.fig(f"struct_{r['name']}", r["name"], "mono", {"text": "\n\n".join(txt) or "No complementary stretch found."},
              how="Top strand 5′→3′; the partner is written 3′→5′ underneath. | marks the reported run of pairs, · other pairs at that offset.",
              yours=" ".join(verdict))
    pair = None
    if len(rows) == 2:
        pair = (rows[0], rows[1])
    elif args.get("pair") and len(args["pair"]) == 2:
        byname = {r["name"]: r for r in rows}
        if all(n in byname for n in args["pair"]):
            pair = (byname[args["pair"][0]], byname[args["pair"][1]])
    if len(rows) > 1:
        cross_art = []
        worst = None
        for x, y in itertools.combinations(rows, 2):
            cd = comp_scan(x["seq"], y["seq"])
            if worst is None or cd["three"]["run"] > worst[2]["three"]["run"]:
                worst = (x, y, cd)
            if cd["three"]["run"] >= 3 or cd["longest"]["run"] >= 8:
                R.flag("warn", f"Cross-dimer <b>{esc(x['name'])}</b> × <b>{esc(y['name'])}</b>: "
                               f"{cd['three']['run']} bp 3′-anchored, {cd['longest']['run']} bp longest run.")
            if len(rows) <= 6 or cd["three"]["run"] >= 3:
                cross_art.append(f"{x['name']} × {y['name']}: longest {cd['longest']['run']} bp, 3′-anchored {cd['three']['run']} bp\n{cd['art']}")
        R.fig("cross", "Cross-dimers", "mono", {"text": "\n\n".join(cross_art) or "No pair shows notable complementarity."},
              how="Every pair of primers slid against each other (antiparallel), as for self-dimers.",
              yours=(f"The strongest 3′-anchored cross-complementarity is {worst[2]['three']['run']} bp "
                     f"({esc(worst[0]['name'])} × {esc(worst[1]['name'])})"
                     + (" — a CHECK: those 3′ ends can prime on each other." if worst[2]["three"]["run"] >= 3 else " — below the 3 bp warning level.")))

    if pair and pair[0]["tm_nn"] is not None and pair[1]["tm_nn"] is not None:
        a, c = pair
        d = abs(a["tm_nn"] - c["tm_nn"])
        low = min(a["tm_nn"], c["tm_nn"])
        R.section("pair", "Pair", f"{esc(a['name'])} + {esc(c['name'])}", "")
        R.table("pair", "Pair analysis", ["", "Value", "Comment"], [
            ["Tm difference", f"{d:.1f} °C", "fine (≤ 5 °C)" if d <= 5 else "CHECK: > 5 °C — the higher-Tm primer will prime less specifically at the lower primer's temperature"],
            ["Suggested annealing (Taq)", f"{low - 5:.0f}–{low - 3:.0f} °C", "lower Tm minus 3–5 °C; confirm with a gradient"],
            ["High-fidelity polymerases", "use the vendor calculator", "Q5/Phusion anneal higher (often Tm + 3 °C with their own Tm values)"],
        ], how="The annealing temperature is a starting point: the lower Tm minus 3–5 °C for Taq-type polymerases in the buffer above.",
            yours=(f"Tm {a['tm_nn']:.1f} and {c['tm_nn']:.1f} °C, {d:.1f} °C apart. Start Taq annealing around {low - 4:.0f} °C."
                   + (" <b>The Tm difference is above 5 °C</b>: consider lengthening the lower-Tm primer at its 5′ end." if d > 5 else "")),
            wide=False)
        if d > 5:
            R.flag("warn", f"Pair {esc(a['name'])} + {esc(c['name'])}: Tm differ by {d:.1f} °C (> 5 °C).")
        else:
            R.flag("ok", f"Pair {esc(a['name'])} + {esc(c['name'])}: Tm within {d:.1f} °C; start Taq annealing near {low - 4:.0f} °C.")

    tpl = args.get("template")
    if tpl:
        from .core_seq import record_from_text, record_json
        if isinstance(tpl, str):
            rec, _ = record_from_text(tpl, bool(args.get("template_circular")), "Pasted template")
            tpl = record_json(rec)
        res = pcr(rows, tpl)
        R.r["pcr"] = {"template": tpl["name"], "sites": res["sites"], "products": res["products"][:50]}
        R.r["record"] = tpl
        prods = res["products"]
        main = [p for p in prods if p["perfect"]]
        R.section("pcr", "Virtual PCR", f"PCR on {esc(tpl['name'])}",
                  f"Binding sites on both strands of the {'circular' if tpl['circular'] else 'linear'} template: at most 2 mismatches "
                  "and the last 3 bases at the 3′ end matching exactly. Primers whose 5′ part does not match (cloning tails) are "
                  "placed by their 3′ 20 nt, and the tail is added to the product size.")
        primer_feats = [{"name": s["primer"], "type": "primer_bind", "strand": s["strand"],
                         "parts": [[s["start0"], s["end0"]]] if s["end0"] > s["start0"] else [[s["start0"], tpl["length"]], [0, s["end0"]]],
                         "start0": s["start0"], "end0": s["end0"], "wraps": s["end0"] <= s["start0"], "length": s["len"],
                         "mm": s["mm"], "qualifiers": {}} for s in res["sites"]]
        R.fig("pcrmap", "Primers on the template", "map", {"record_ref": True, "primers": primer_feats,
                                                          "products": [p for p in prods[:6]], "cutters": []},
              how="Arrows are primer binding sites (pointing in the direction of extension); the shaded arc is the predicted product.",
              yours=(f"{pl(len(res['sites']), 'binding site')} for {pl(len(rows), 'primer')}; {pl(len(prods), 'product')} predicted"
                     + (f", the main one {commas(main[0]['size'])} bp from {esc(main[0]['fwd'])} to {esc(main[0]['rev'])}"
                        + (" (across the origin)" if main[0]["wraps"] else "") if main else "") + "."),
              wide=True)
        R.table("products", "Predicted products", ["Forward", "Reverse", "Size (bp)", "Template span (1-based)", "Mismatches", "Note"],
                [[p["fwd"], p["rev"], commas(p["size"]), f"{commas(p['start0'] + 1)}–{commas(p['end0'])}" + (" (across origin)" if p["wraps"] else ""),
                  p["mm"], "main product" if p is (main[0] if main else None) else ("off-target" if p["mm"] else "alternative product")]
                 for p in prods[:50]],
                how="Every forward site paired with every downstream reverse site within 15 kb (on a circle, products can run through the origin).",
                yours=("<b>No product</b>: no forward/reverse pair of sites faces each other on this template. Check that the primers "
                       "are written 5′→3′ and that one is the reverse complement of the top strand." if not prods else
                       f"{pl(len(prods), 'product')}; " + (f"{pl(len(prods) - 1, 'additional product')} — check for extra bands."
                                                           if len(prods) > 1 else "a single product, as intended.")))
        R.table("sites", "Binding sites (on-target and off-target)", ["Primer", "Strand", "Start–end (1-based)", "Mismatches", "5′ tail (nt)", "Inside"],
                [[s["primer"], "+" if s["strand"] == 1 else "−", f"{commas(s['start0'] + 1)}–{commas(s['end0'])}", s["mm"], s["tail"],
                  ", ".join(s["features"]) or "—"] for s in res["sites"]],
                how="A site with 1–2 mismatches (3′ end perfect) can still prime, particularly at low annealing temperatures.",
                yours=_sites_yours(rows, res["sites"]))
        unbound = [r["name"] for r in rows if not any(s["primer"] == r["name"] for s in res["sites"])]
        if unbound:
            R.flag("warn", f"{', '.join(esc(u) for u in unbound)} do{'es' if len(unbound) == 1 else ''} not bind the template "
                           "(≤ 2 mismatches with a perfect 3′ end), on either strand.")
        if prods and not main:
            R.flag("warn", "Every predicted product needs a mismatched primer; there is no perfect-match product.")
        elif main:
            R.flag("ok", f"Main product {commas(main[0]['size'])} bp ({esc(main[0]['fwd'])} → {esc(main[0]['rev'])}).")
        offs = [s for s in res["sites"] if s["mm"] > 0]
        if offs:
            R.flag("info", f"{pl(len(offs), 'off-target site')} with 1–2 mismatches — listed in the binding-site table.")
    R.method("Melting temperature", "Biopython <code>Bio.SeqUtils.MeltingTemp</code>: <code>Tm_NN(seq, nn_table=DNA_NN3, "
             f"saltcorr=7, Na={b['Na']:g}, K={b['K']:g}, Tris={b['Tris']:g}, Mg={b['Mg']:g}, dNTPs={b['dNTPs']:g}, dnac1={b['primer_nM']:g}, "
             f"dnac2=0)</code>" + (f" followed by <code>chem_correction(DMSO={b['DMSO']:g})</code>" if b["DMSO"] else "")
             + ". Nearest-neighbour parameters: Allawi & SantaLucia, Biochemistry 1997; SantaLucia, PNAS 1998. Salt correction: "
               "Owczarzy et al., Biochemistry 2008. Also <code>Tm_GC(valueset=7)</code> and <code>Tm_Wallace</code> (Wallace et al. 1979).")
    R.method("Molecular weight", "<code>Bio.SeqUtils.molecular_weight(seq, 'DNA')</code>, single-stranded, with 5′ and 3′ hydroxyls "
                                 "(no 5′ phosphate, no dye).")
    R.method("Complementarity", "A base-pairing heuristic written for this app (longest Watson–Crick run between antiparallel "
                                "strands, 3′-anchored runs, hairpin stems ≥ 4 bp with loops ≥ 3 nt). No free energies are computed.")
    if tpl:
        R.method("Virtual PCR", "Exhaustive window scan of both strands; ≤ 2 mismatches, the 3′-terminal 3 bases exact; products "
                                "≤ 15 kb between facing sites.")
    return R.done()


def _sites_yours(rows, sites):
    parts = []
    for r in rows:
        mine = [s for s in sites if s["primer"] == r["name"]]
        perfect = [s for s in mine if s["mm"] == 0]
        parts.append(f"{esc(r['name'])}: {pl(len(perfect), 'perfect site')}"
                     + (f", {len(mine) - len(perfect)} with mismatches" if len(mine) > len(perfect) else ""))
    return "; ".join(parts) + "."
