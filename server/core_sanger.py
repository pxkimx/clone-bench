"""Sanger verification: read .ab1 traces, Mott-trim, align to a reference on both strands, call every difference
with a confidence, find mixed peaks, and say which parts of the construct are verified. Biopython + stdlib."""
from __future__ import annotations

import io
import re

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.Seq import Seq

from .core import Result, UserFacingError, commas, esc, file_bytes, pl
from .core_seq import (check_cds, extract_feature, fmt_span, record_from_text, record_json, parse_sequence_file,
                       revcomp)

CHANNEL_KEYS = ["DATA9", "DATA10", "DATA11", "DATA12"]
RAW_KEYS = ["DATA1", "DATA2", "DATA3", "DATA4"]
MIXED_RATIO = 0.30


# ---------------------------------------------------------------- reading a trace
def read_trace(f: dict) -> dict:
    name = f.get("name") or "read.ab1"
    data = file_bytes(f)
    if not data:
        raise UserFacingError(f"{name} is empty.")
    if data[:4] != b"ABIF":
        low = name.lower()
        if data[:1] in (b">", b"L") or low.endswith((".gb", ".gbk", ".fa", ".fasta", ".dna")):
            raise UserFacingError(f"{name} is a sequence file, not a Sanger trace. Load it as the reference, and the .ab1 "
                                  "file(s) as reads.")
        if data[:4] == b".scf":
            raise UserFacingError(f"{name} is an SCF trace; this app reads ABI .ab1 files. Most sequencing providers can "
                                  "send .ab1 instead.")
        raise UserFacingError(f"{name} is not an ABI trace file (.ab1): it does not start with the ABIF signature.")
    try:
        rec = SeqIO.read(io.BytesIO(data), "abi")
    except Exception as e:  # noqa: BLE001
        raise UserFacingError(f"{name} could not be read as an ABI trace: {str(e)[:200]}. The file may be truncated.")
    raw = rec.annotations.get("abif_raw", {})
    order = raw.get("FWO_1", b"GATC")
    order = order.decode() if isinstance(order, bytes) else str(order)
    notes = []
    keys = CHANNEL_KEYS if all(k in raw for k in CHANNEL_KEYS) else None
    if keys is None and all(k in raw for k in RAW_KEYS):
        keys = RAW_KEYS
        notes.append(("info", f"{esc(name)} has no processed trace channels (DATA9–12); the raw channels DATA1–4 are drawn "
                              "instead, so peaks look noisier and less evenly spaced than the basecaller saw them."))
    channels = {}
    if keys:
        for k, base in zip(keys, order):
            channels[base] = [int(v) for v in raw[k]]
    else:
        notes.append(("warn", f"{esc(name)} contains no trace channels at all — only the basecalls can be analysed."))
    ploc = list(raw.get("PLOC2") or raw.get("PLOC1") or [])
    calls = str(rec.seq).upper()
    quals = list(rec.letter_annotations.get("phred_quality") or [])
    if len(quals) != len(calls):
        quals = [0] * len(calls)
    has_q = any(q > 0 for q in quals)
    if len(ploc) != len(calls):
        ploc = ploc[:len(calls)] if len(ploc) > len(calls) else []
    if not calls:
        raise UserFacingError(f"{name} contains no basecalls.")
    sample = raw.get("SMPL1") or rec.id
    sample = sample.decode("latin-1", "replace") if isinstance(sample, bytes) else str(sample)
    return {"name": name, "sample": sample, "calls": calls, "quals": quals, "has_q": has_q, "ploc": ploc,
            "channels": channels, "order": order, "notes": notes, "raw_channels": keys == RAW_KEYS}


def mott_trim(quals: list, cutoff: float = 0.05, segment: int = 20) -> tuple[int, int]:
    """Richard Mott's modified trimming, the algorithm Biopython's 'abi-trim' format uses (Bio.SeqIO.AbiIO._abi_trim):
    each base scores cutoff − 10^(−Q/10); keep the segment with the highest running sum. Returns [start, end)."""
    if len(quals) <= segment:
        return 0, len(quals)
    score = [cutoff - 10 ** (q / -10.0) for q in quals]
    cum = [0.0]
    start, started = 0, False
    for i in range(1, len(score)):
        v = cum[-1] + score[i]
        if v < 0:
            cum.append(0.0)
        else:
            cum.append(v)
            if not started:
                start, started = i, True
    end = cum.index(max(cum))
    return start, end


# ---------------------------------------------------------------- alignment
def aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -6
    a.extend_gap_score = -1
    return a


def align_read(read: str, ref: str, circular: bool) -> dict | None:
    """Local alignment of the (trimmed) read to the reference on both strands; circular references are
    aligned as reference + reference so a read across the origin aligns in one piece."""
    tgt = ref + ref if circular else ref
    al = aligner()
    fwd = al.score(tgt, read)
    rev = al.score(tgt, revcomp(read))
    strand = 1 if fwd >= rev else -1
    q = read if strand == 1 else revcomp(read)
    best = max(fwd, rev)
    if best < 40:          # fewer than ~20 matching bases: not this reference
        return None
    aln = al.align(tgt, q)[0]
    return {"strand": strand, "q": q, "score": float(best), "coords": [[int(x) for x in row] for row in aln.coordinates],
            "other_score": float(min(fwd, rev))}


def run_len(s: str, i: int, step: int) -> int:
    """Length of the run of identical bases starting at i and going in direction step."""
    if not 0 <= i < len(s):
        return 0
    b, n = s[i], 0
    while 0 <= i < len(s) and s[i] == b:
        n += 1
        i += step
    return n


def homopolymer_at(ref: str, p: int) -> int:
    """Longest run of one base touching position p (either side)."""
    L = len(ref)
    best = 0
    for j in (p - 1, p, p + 1):
        if 0 <= j < L:
            best = max(best, run_len(ref, j, -1) + run_len(ref, j, 1) - 1)
    return best


# ---------------------------------------------------------------- CDS consequences
AA3 = dict(protein_letters_1to3)
AA3["*"] = "Ter"
AA3["X"] = "Xaa"


def cds_index(fj: dict, p: int):
    """0-based index into the CDS (as read 5′→3′ on its strand) of reference base p, or None."""
    acc = 0
    for s, e in fj["parts"]:
        if s <= p < e:
            return acc + ((p - s) if fj["strand"] != -1 else (e - 1 - p))
        acc += e - s
    return None


def codon_change(fj, ref, p, alt_base, table=11):
    """Consequence of replacing reference base p with alt_base inside CDS feature fj."""
    k = cds_index(fj, p)
    if k is None:
        return None
    cs = int((fj.get("qualifiers") or {}).get("codon_start") or 1) - 1
    k -= cs
    if k < 0:
        return None
    try:
        table = int((fj.get("qualifiers") or {}).get("transl_table") or table)
    except ValueError:
        table = 11
    cds = extract_feature(ref, fj)[cs:]
    b = alt_base if fj["strand"] != -1 else revcomp(alt_base)
    if b not in "ACGT":
        return {"feature": fj["name"], "codon": k // 3 + 1, "kind": "ambiguous", "text": "ambiguous basecall — consequence not assessed"}
    i = k // 3
    rc = cds[3 * i:3 * i + 3]
    if len(rc) < 3:
        return {"feature": fj["name"], "codon": i + 1, "kind": "partial", "text": "in an incomplete final codon"}
    mc = rc[:k % 3] + b + rc[k % 3 + 1:]
    ra = str(Seq(rc).translate(table=table))
    ma = str(Seq(mc).translate(table=table))
    from Bio.Data import CodonTable
    starts = CodonTable.unambiguous_dna_by_id[table].start_codons
    if i == 0 and rc in starts:
        ra = "M"
        if mc in starts:
            ma = "M"
    n = i + 1
    if i == 0 and ra == "M" and ma != "M":
        kind, txt = "start-loss", f"start codon lost (p.Met1?; {rc}→{mc})"
    elif ra == ma:
        kind, txt = "synonymous", f"synonymous (p.{AA3[ra]}{n}=, {ra}{n}=; {rc}→{mc})"
    elif ma == "*":
        kind, txt = "nonsense", f"nonsense (p.{AA3[ra]}{n}Ter, {ra}{n}*; {rc}→{mc})"
    elif ra == "*":
        kind, txt = "stop-loss", f"stop codon lost (p.Ter{n}{AA3[ma]}ext; {rc}→{mc})"
    else:
        kind, txt = "missense", f"missense p.{AA3[ra]}{n}{AA3[ma]} ({ra}{n}{ma}; {rc}→{mc})"
    return {"feature": fj["name"], "codon": n, "kind": kind, "text": txt, "short": f"{ra}{n}{ma}" if kind == "missense" else None,
            "hgvs": f"p.{AA3[ra]}{n}{AA3[ma]}" if kind == "missense" else None}


def indel_effect(fj, p, n_bases, kind):
    k = cds_index(fj, p if kind == "deletion" else max(p - 1, 0))
    if k is None:
        return None
    codon = k // 3 + 1
    if n_bases % 3:
        return {"feature": fj["name"], "codon": codon, "kind": "frameshift",
                "text": f"frameshift ({n_bases}-bp {kind} at codon {codon}; everything downstream is out of frame)"}
    return {"feature": fj["name"], "codon": codon, "kind": "in-frame",
            "text": f"in-frame {kind} of {n_bases // 3} codon{'s' if n_bases > 3 else ''} at codon {codon}"}


# ---------------------------------------------------------------- one read
def analyse_read(tr: dict, rj: dict) -> dict:
    ref, L, circ = rj["seq"], rj["length"], rj["circular"]
    calls, quals, has_q = tr["calls"], tr["quals"], tr["has_q"]
    ts, te = mott_trim(quals) if has_q else (0, len(calls))
    trim_note = None
    if has_q and te - ts < 30:
        trim_note = "Mott trimming left fewer than 30 good bases; the whole read was aligned instead."
        ts, te = 0, len(calls)
    read = calls[ts:te]
    q_all = [q >= 20 for q in quals]
    stats = {"length": len(calls), "trim_start": ts, "trim_end": te, "trimmed_length": te - ts,
             "q20": 100 * sum(q_all) / len(calls) if has_q else None,
             "q30": 100 * sum(q >= 30 for q in quals) / len(calls) if has_q else None,
             "mean_q_trimmed": sum(quals[ts:te]) / max(1, te - ts) if has_q else None,
             "q20_trimmed": 100 * sum(q >= 20 for q in quals[ts:te]) / max(1, te - ts) if has_q else None}
    out = {"name": tr["name"], "sample": tr["sample"], "stats": stats, "has_q": has_q, "trim_note": trim_note,
           "diffs": [], "mixed": [], "aligned": False}
    aln = align_read(read, ref, circ)
    if aln is None:
        return out
    strand, q, coords = aln["strand"], aln["q"], aln["coords"]
    n = te - ts

    def orig(qi):           # index in q -> index in the original read
        return ts + (qi if strand == 1 else n - 1 - qi)

    def Q(i):
        return quals[i] if 0 <= i < len(quals) else 0

    def minQ(i, r=2):
        vals = [Q(j) for j in range(i - r, i + r + 1) if 0 <= j < len(quals) and j != i]
        return min(vals) if vals else 0

    cov = {}               # ref position -> state (1 low-Q, 2 verified, 3 confident difference)
    qmap = {}              # q index -> ref position (aligned bases only)
    diffs = []
    ident = mism = 0
    t_start, t_end = coords[0][0], coords[0][-1]
    for k in range(len(coords[0]) - 1):
        t0, t1 = coords[0][k], coords[0][k + 1]
        q0, q1 = coords[1][k], coords[1][k + 1]
        if t1 > t0 and q1 > q0:
            for d in range(t1 - t0):
                tp, qi = t0 + d, q0 + d
                rp = tp % L
                oi = orig(qi)
                qmap[qi] = rp
                rb, cb = ref[rp], q[qi]
                good = Q(oi) >= 20 or not has_q
                if rb == cb:
                    ident += 1
                    cov[rp] = max(cov.get(rp, 0), 2 if (good and has_q) else 1)
                else:
                    mism += 1
                    diffs.append({"type": "substitution", "ref_pos": rp, "ref": rb, "read": cb, "read_index": oi,
                                  "q": Q(oi), "q_min": minQ(oi), "q_idx": qi})
                    cov[rp] = max(cov.get(rp, 0), 3 if Q(oi) >= 20 else 1)
        elif t1 > t0:      # reference bases missing from the read: a deletion in the read
            left, right = orig(q0 - 1) if q0 > 0 else None, orig(q0) if q0 < len(q) else None
            qs = [Q(i) for i in (left, right) if i is not None]
            qv = min(qs) if qs else 0
            rp = t0 % L
            diffs.append({"type": "deletion", "ref_pos": rp, "ref": "".join(ref[(t0 + d) % L] for d in range(t1 - t0)), "read": "-",
                          "read_index": min(i for i in (left, right) if i is not None) if qs else None, "q": qv,
                          "q_min": min(minQ(i) for i in (left, right) if i is not None) if qs else 0, "len": t1 - t0, "q_idx": q0})
            for d in range(t1 - t0):
                cov[(t0 + d) % L] = max(cov.get((t0 + d) % L, 0), 3 if qv >= 20 else 1)
        elif q1 > q0:      # read bases absent from the reference: an insertion in the read
            idx = [orig(qi) for qi in range(q0, q1)]
            qv = min(Q(i) for i in idx)
            rp = t0 % L
            diffs.append({"type": "insertion", "ref_pos": rp, "ref": "-", "read": q[q0:q1], "read_index": min(idx), "q": qv,
                          "q_min": min(minQ(i) for i in idx), "len": q1 - q0, "q_idx": q0})
    # annotate each difference
    for d in diffs:
        p = d["ref_pos"]
        hp = homopolymer_at(ref, p) if d["type"] != "substitution" else max(run_len(ref, p - 1, -1), run_len(ref, (p + 1) % L, 1))
        if d["type"] != "substitution":
            base = (d["ref"] if d["type"] == "deletion" else d["read"])
            hp_indel = len(set(base)) == 1 and (ref[p % L] == base[0] or ref[(p - 1) % L] == base[0])
            d["homopolymer"] = hp if hp_indel and hp >= 3 else 0
        else:
            d["homopolymer"] = hp if hp >= 4 else 0
        ambiguous = d["type"] == "substitution" and d["read"] not in "ACGT"
        if not has_q:
            conf = "unscored"
        elif d["q"] < 20:
            conf = "probably a basecall error"
        elif ambiguous or d["q"] < 30 or d["q_min"] < 20 or d["homopolymer"]:
            conf = "check the trace"
        else:
            conf = "confident"
        d["confidence"] = conf
        d["ambiguous"] = ambiguous
        feats, effects = [], []
        for fj in rj["features"]:
            inside = any(s <= p < e for s, e in fj["parts"])
            if not inside:
                continue
            feats.append(fj["name"])
            if fj["type"] == "CDS":
                if d["type"] == "substitution":
                    eff = codon_change(fj, ref, p, d["read"])
                else:
                    eff = indel_effect(fj, p, d["len"], d["type"])
                if eff:
                    effects.append(eff)
        d["features"] = feats
        d["effects"] = effects
        d["ref_pos1"] = p + 1
        d["read_pos1"] = (d["read_index"] + 1) if d.get("read_index") is not None else None
        d["trace_x"] = tr["ploc"][d["read_index"]] if tr["ploc"] and d.get("read_index") is not None and d["read_index"] < len(tr["ploc"]) else None
        if d["type"] == "deletion" and d["trace_x"] is not None and d["read_index"] + 1 < len(tr["ploc"]):
            # a missing base sits between two called peaks: mark the gap, not the peak before it
            d["trace_x"] = (tr["ploc"][d["read_index"]] + tr["ploc"][d["read_index"] + 1]) / 2
    out["diffs"] = diffs
    covered_ref = (t_end - t_start)
    out.update({"aligned": True, "strand": strand, "score": aln["score"], "ref_start": t_start % L, "ref_end": (t_end - 1) % L + 1,
                "ref_span": covered_ref, "wraps": circ and t_start % L > (t_end - 1) % L,
                "identity": 100 * ident / max(1, ident + mism + sum(d.get("len", 1) for d in diffs if d["type"] != "substitution")),
                "coverage": _runs(cov, L)})
    out["alignment_text"] = alignment_rows(ref, q, coords, L, strand, orig, quals, has_q)
    # mixed peaks
    out["mixed"] = mixed_peaks(tr, ts, te, qmap, strand, n)
    return out


def _runs(cov: dict, L: int) -> list:
    """Run-length [start0, end0, state] over the reference, state 0 = not covered."""
    runs = []
    cur, s0 = None, 0
    for p in range(L):
        st = cov.get(p, 0)
        if st != cur:
            if cur is not None:
                runs.append([s0, p, cur])
            cur, s0 = st, p
    runs.append([s0, L, cur])
    return runs


def alignment_rows(ref, q, coords, L, strand, orig, quals, has_q, width=60):
    """Blocks of three lines (reference, match bars, read) with positions, for the monospace view."""
    top, mid, bot, refpos, readidx = [], [], [], [], []
    for k in range(len(coords[0]) - 1):
        t0, t1 = coords[0][k], coords[0][k + 1]
        q0, q1 = coords[1][k], coords[1][k + 1]
        if t1 > t0 and q1 > q0:
            for d in range(t1 - t0):
                rb, cb = ref[(t0 + d) % L], q[q0 + d]
                top.append(rb); bot.append(cb); mid.append("|" if rb == cb else "x")
                refpos.append((t0 + d) % L); readidx.append(orig(q0 + d))
        elif t1 > t0:
            for d in range(t1 - t0):
                top.append(ref[(t0 + d) % L]); bot.append("-"); mid.append(" ")
                refpos.append((t0 + d) % L); readidx.append(None)
        elif q1 > q0:
            for d in range(q1 - q0):
                top.append("-"); bot.append(q[q0 + d]); mid.append(" ")
                refpos.append(None); readidx.append(orig(q0 + d))
    return {"ref": "".join(top), "mid": "".join(mid), "read": "".join(bot), "refpos": refpos, "readidx": readidx,
            "strand": strand}


def mixed_peaks(tr, ts, te, qmap, strand, n, ratio=MIXED_RATIO):
    ch, ploc, calls, quals = tr["channels"], tr["ploc"], tr["calls"], tr["quals"]
    if not ch or not ploc or len(ch) < 4:
        return []
    bases = list(ch)
    heights = []
    for i in range(ts, te):
        x = ploc[i]
        vals = {b: max(ch[b][max(0, x - 2):x + 3] or [0]) for b in bases}
        heights.append(vals)
    prim_all = sorted(max(v.values()) for v in heights)
    floor = prim_all[len(prim_all) // 2] * 0.15 if prim_all else 0
    out = []
    for i, vals in zip(range(ts, te), heights):
        call = calls[i]
        prim_b = call if call in vals else max(vals, key=vals.get)
        prim = vals[prim_b]
        if prim <= 0 or prim < floor:
            continue
        sec_b = max((b for b in bases if b != prim_b), key=lambda b: vals[b])
        r = vals[sec_b] / prim
        if r >= ratio:
            qi = (i - ts) if strand == 1 else (n - 1 - (i - ts))
            rp = qmap.get(qi)
            out.append({"read_index": i, "read_pos1": i + 1, "call": call, "primary": prim_b, "secondary": sec_b,
                        "ratio": round(r, 2), "q": quals[i], "ref_pos1": rp + 1 if rp is not None else None,
                        "trace_x": ploc[i]})
    return out


# ---------------------------------------------------------------- reference & demo
def reference_from(args) -> dict:
    ref = args.get("reference")
    if isinstance(ref, dict) and ref.get("seq"):
        return ref
    if args.get("reference_file"):
        rec, _, _ = parse_sequence_file(args["reference_file"], None if args.get("reference_circular") in (None, "", "auto")
                                        else bool(args["reference_circular"]))
        rj = record_json(rec)
    elif args.get("reference_text"):
        rec, _ = record_from_text(args["reference_text"], bool(args.get("reference_circular")), "Pasted reference")
        rj = record_json(rec)
    else:
        raise UserFacingError("Add a reference: pick a loaded construct, paste a sequence, or drop a GenBank/FASTA file.")
    return rj


def build_demo_reference(ab1: dict) -> dict:
    """The honest demo: a reference made from 3730.ab1's own Mott-trimmed basecalls with two planted edits — one
    substitution in a high-quality stretch and one extra base (so the read shows a 1-bp deletion) — plus the ORF the
    read contains annotated as a CDS, so the consequence calls can be seen."""
    tr = read_trace(ab1)
    ts, te = mott_trim(tr["quals"])
    calls, q = tr["calls"][ts:te], tr["quals"][ts:te]
    ref = list(calls)
    # replace ambiguous calls in the reference by the strongest channel's base, so the reference is plain DNA
    for i, c in enumerate(ref):
        if c not in "ACGT":
            x = tr["ploc"][ts + i]
            ref[i] = max(tr["channels"], key=lambda b: tr["channels"][b][x])
    ref = "".join(ref)
    atg = ref.find("ATG")
    atg = ref.find("ATGGTGAGCAAGG") if ref.find("ATGGTGAGCAAGG") >= 0 else atg
    good = lambda i, r=4: all(q[j] >= 40 for j in range(i - r, i + r + 1))  # noqa: E731

    def missense_at(i):
        k = i - atg
        cod = ref[atg + k // 3 * 3: atg + k // 3 * 3 + 3]
        alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref[i]]
        mut = cod[:k % 3] + alt + cod[k % 3 + 1:]
        a, b = str(Seq(cod).translate()), str(Seq(mut).translate())
        return alt if a != b and "*" not in a + b else None

    sub = next(i for i in range(atg + 150, len(ref) - 300) if (i - atg) % 3 == 1 and good(i) and missense_at(i))
    alt = missense_at(sub)
    # the insertion site: high quality, and the inserted base differs from both neighbours (no homopolymer ambiguity)
    ins = next(i for i in range(sub + 240, len(ref) - 100) if good(i) and ref[i - 1] != ref[i])
    extra = next(b for b in "ACGT" if b not in (ref[i - 1] for i in [ins]) and b != ref[ins])
    planted = ref[:sub] + alt + ref[sub + 1:ins] + extra + ref[ins:]
    orf_end = atg + (len(planted) - atg) // 3 * 3
    return {"seq": planted, "atg": atg, "orf_end": orf_end, "substitution": {"ref_pos1": sub + 1, "ref": alt, "read": ref[sub]},
            "extra_base": {"ref_pos1": ins + 1, "base": extra}, "trim": [ts, te]}


def demo_reference_api(args):
    d = build_demo_reference(args["file"])
    return d


def demo_record(d: dict) -> dict:
    """GenBank-like record JSON for the demo reference."""
    L = len(d["seq"])
    return {"name": "demo_reference", "id": "demo_reference", "length": L, "circular": False, "seq": d["seq"],
            "description": "3730.ab1 basecalls with two planted edits",
            "features": [{"i": 0, "name": "ORF (mCherry-linker-GFP start, as read)", "type": "CDS", "strand": 1,
                          "parts": [[d["atg"], d["orf_end"]]], "start0": d["atg"], "end0": d["orf_end"], "wraps": False,
                          "length": d["orf_end"] - d["atg"], "fuzzy": True, "partial5": False, "partial3": True,
                          "qualifiers": {"note": "ORF in this read's own basecalls; runs off the end of the read",
                                         "transl_table": "11"}}]}


# ---------------------------------------------------------------- the tool
def sanger(args: dict) -> dict:
    files = args.get("files") or []
    if not files:
        raise UserFacingError("Drop one or more .ab1 trace files.")
    if len(files) > 24:
        raise UserFacingError("Load at most 24 reads at a time.")
    rj = reference_from(args)
    ref, L, circ = rj["seq"], rj["length"], rj["circular"]
    for f in rj["features"]:
        if f["type"] == "CDS" and "translation" not in f:
            f["translation"] = check_cds(f, ref)["translation"]
    traces = [read_trace(f) for f in files]
    reads = [analyse_read(t, rj) for t in traces]
    focus = None
    if args.get("focus_feature") not in (None, ""):
        from .core_seq import find_feature
        focus = find_feature(rj, args["focus_feature"])
    name = (traces[0]["name"] if len(traces) == 1 else f"{len(traces)} reads") + f" vs {rj['name']}"
    R = Result("sanger", name)
    R.r["record"] = rj
    R.r["demo"] = bool(args.get("demo"))
    if args.get("demo"):
        R.flag("info", "<b>Demo</b>: the reference was built from this read's own high-quality basecalls with two planted edits "
                       "(one substitution, one extra base), so you can see what a difference looks like. It is not a real clone.")
    for t in traces:
        for lvl, txt in t["notes"]:
            R.flag(lvl, txt)
    aligned = [r for r in reads if r["aligned"]]
    for r in reads:
        if not r["aligned"]:
            R.flag("warn", f"<b>{esc(r['name'])}</b> does not align to {esc(rj['name'])} on either strand. Is it the right "
                           "reference, or a failed read (check the trace for signal)?")
        if not r["has_q"]:
            R.flag("warn", f"<b>{esc(r['name'])}</b> carries no quality values (the basecaller wrote none — common for old 310 "
                           "files). Nothing was trimmed and no difference can be called confident: every one needs a look at the trace.")
        if r.get("trim_note"):
            R.flag("warn", f"<b>{esc(r['name'])}</b>: {r['trim_note']}")
    conf = [d for r in aligned for d in r["diffs"] if d["confidence"] == "confident"]
    check = [d for r in aligned for d in r["diffs"] if d["confidence"] in ("check the trace", "unscored")]
    R.tile(len(reads), "reads")
    if aligned:
        R.tile(f"{min(r['identity'] for r in aligned):.1f}%" if len(aligned) > 1 else f"{aligned[0]['identity']:.1f}%", "identity (aligned region)")
    R.tile(len(conf), "confident differences")
    R.tile(len(check), "to check in the trace")

    # combined coverage
    comb = [0] * L
    for r in aligned:
        for s, e, st in r["coverage"]:
            for p in range(s, e):
                if st > comb[p] or (st == 3 and comb[p] != 3):
                    comb[p] = max(comb[p], st)
    comb_runs = _runs({p: v for p, v in enumerate(comb) if v}, L)
    verified = sum(1 for v in comb if v == 2)
    R.tile(f"{100 * verified / L:.1f}%", "of reference verified at Q ≥ 20")
    R.r["summary"] = {"reads": len(reads), "confident": len(conf), "verified_pct": 100 * verified / L}

    # per-read summaries (the "yours" sentences)
    for r in aligned:
        rc = [d for d in r["diffs"] if d["confidence"] == "confident"]
        rk = [d for d in r["diffs"] if d["confidence"] in ("check the trace", "unscored")]
        rn = [d for d in r["diffs"] if d["confidence"] == "probably a basecall error"]
        span = f"bp {commas(r['ref_start'] + 1)}–{commas(r['ref_end'])}" + (" (across the origin)" if r["wraps"] else "")
        s = f"<b>{esc(r['name'])}</b> covers {span} of {esc(rj['name'])} at {r['identity']:.1f} % identity ({'forward' if r['strand'] == 1 else 'reverse'} strand). "
        if not r["diffs"]:
            s += "No differences from the reference."
        else:
            if rc:
                s += f"{pl(len(rc), 'confident difference')}: " + "; ".join(diff_phrase(d) for d in rc[:4]) + ". "
            if rk:
                s += f"{pl(len(rk), 'difference')} to check in the trace. "
            if rn:
                s += f"{pl(len(rn), 'difference')} at low quality (Q < 20) {'is' if len(rn) == 1 else 'are'} probably not real."
        r["yours"] = s
        R.flag("warn" if rc else "ok", s)

    # ---- per read sections
    for k, (t, r) in enumerate(zip(traces, reads)):
        sid = f"read{k}"
        R.section(sid, f"Read {k + 1}", esc(t["name"]), f"Sample {esc(t['sample'])} · {commas(len(t['calls']))} bases called"
                  + (f" · trimmed to {commas(r['stats']['trimmed_length'])} (Mott, cutoff 0.05)" if r["has_q"] else " · no quality values"))
        st = r["stats"]
        R.fig(f"trace{k}", "Chromatogram", "trace",
              {"read": k, "channels": t["channels"], "order": t["order"], "ploc": t["ploc"], "calls": t["calls"], "quals": t["quals"],
               "trim": [st["trim_start"], st["trim_end"]], "diffs": [{"x": d["trace_x"], "i": d.get("read_index"), "type": d["type"],
                                                                     "confidence": d["confidence"], "label": diff_label(d)} for d in r["diffs"]],
               "mixed": [{"x": m["trace_x"], "i": m["read_index"], "ratio": m["ratio"]} for m in r["mixed"]],
               "raw": t["raw_channels"]},
              how=("The four fluorescence channels from the .ab1 file (processed DATA9–12, in the dye order given by FWO_1): "
                   "A green, C blue, G black, T red. Letters are the basecalls at the peak positions (PLOC2), bars under them "
                   "the Phred quality (Q20 = 1 % error, Q30 = 0.1 %). Grey shading is what Mott trimming removed. Scroll or drag "
                   "sideways, zoom with the buttons, and jump between differences."),
              yours=(f"{commas(st['length'])} bases; " + (f"{st['q20']:.1f} % at Q ≥ 20 and {st['q30']:.1f} % at Q ≥ 30 over the whole read; "
                                                             f"after trimming {commas(st['trimmed_length'])} bases remain with mean Q {st['mean_q_trimmed']:.1f}."
                                                             if r["has_q"] else "the file has no quality values.")
                     + (f" {pl(len(r['mixed']), 'possible mixed peak')}." if r["mixed"] else "")),
              wide=True, sub="scroll · zoom · jump to difference")
        if r["has_q"]:
            R.fig(f"qual{k}", "Quality along the read", "qual", {"quals": t["quals"], "trim": [st["trim_start"], st["trim_end"]]},
                  how="Phred quality per base with the Q20 and Q30 lines; the shaded ends are what Mott trimming (cutoff 0.05) removed.",
                  yours=(f"Trimmed region {commas(st['trim_start'] + 1)}–{commas(st['trim_end'])}: mean Q {st['mean_q_trimmed']:.1f}, "
                         f"{st['q20_trimmed']:.1f} % of its bases at Q ≥ 20."))
        if not r["aligned"]:
            continue
        R.fig(f"aln{k}", "Alignment to the reference", "alignment", {"read": k, **r["alignment_text"]},
              how=("Local alignment (<code>PairwiseAligner</code>, match 2, mismatch −3, gap open −6, extend −1) of the trimmed read "
                   f"to the reference{' (as reference + reference, so reads across the origin align in one piece)' if circ else ''}, "
                   "on whichever strand scores higher. Differences are highlighted; click one to jump to it in the trace."),
              yours=r["yours"], wide=True)
        rows = []
        for d in r["diffs"]:
            rows.append([d["type"], commas(d["ref_pos1"]), d["read_pos1"] or "—", f"{d['ref']} → {d['read']}",
                         d["q"] if r["has_q"] else "—", d["q_min"] if r["has_q"] else "—", ", ".join(d["features"]) or "—",
                         "; ".join(e["text"] for e in d["effects"]) or "—",
                         d["confidence"] + (f" (homopolymer ×{d['homopolymer']})" if d["homopolymer"] else "")])
        R.table(f"diffs{k}", "Differences", ["Type", "Ref pos (1-based)", "Read pos", "Ref → read", "Q", "Min Q ±2",
                                             "Feature", "Consequence", "Call"], rows,
                how=("One row per substitution, insertion or deletion. <b>confident</b>: Q ≥ 30 with neighbours ≥ 20; <b>check the "
                     "trace</b>: Q 20–30, a weak neighbour, an ambiguous basecall, or next to a homopolymer (indels in runs of one "
                     "base are a classic basecalling error); <b>probably a basecall error</b>: Q < 20. Inside a CDS the codon change "
                     "and its consequence are given in 3-letter (HGVS) and 1-letter form."),
                yours=r["yours"])
        if r["mixed"]:
            R.table(f"mixed{k}", "Possible mixed peaks", ["Read pos", "Ref pos", "Called", "Second peak", "Second / first", "Q"],
                    [[m["read_pos1"], commas(m["ref_pos1"]) if m["ref_pos1"] else "—", m["primary"], m["secondary"], f"{m['ratio']:.2f}", m["q"]]
                     for m in sorted(r["mixed"], key=lambda m: -m["ratio"])[:200]],
                    how=(f"At each called peak the tallest other channel is compared with the called one (within ±2 scans). A second "
                         f"peak above {int(MIXED_RATIO * 100)} % of the first can mean a mixed population (two plasmids in the colony) "
                         "or a heterozygous site — or just noise, which is why they are reported separately from the differences."),
                    yours=mixed_yours(r))

    # ---- coverage
    R.section("coverage", "Verification", "Is the clone fully sequence-verified?",
              "Every reference base is <b>verified</b> when a read covers it with a matching base at Q ≥ 20.")
    R.fig("coverage", "Coverage of the reference", "coverage",
          {"length": L, "circular": circ, "reads": [{"name": r["name"], "runs": r.get("coverage") or [[0, L, 0]], "strand": r.get("strand")} for r in reads],
           "combined": comb_runs, "features": [{"name": f["name"], "type": f["type"], "parts": f["parts"], "strand": f["strand"]} for f in rj["features"]],
           "diffs": [{"pos": d["ref_pos"], "confidence": d["confidence"]} for r in aligned for d in r["diffs"]]},
          how="One track per read, then all reads combined: teal = verified (match at Q ≥ 20), amber = covered at low quality, "
              "red = a difference, grey = not covered.",
          yours=coverage_yours(rj, comb, focus), wide=True)
    frows = []
    for f in rj["features"]:
        pos = [p for s, e in f["parts"] for p in range(s, e)]
        if not pos:
            continue
        v = sum(1 for p in pos if comb[p] == 2)
        gaps = _gaps(pos, comb)
        frows.append([f["name"], f["type"], fmt_span(f), f"{100 * v / len(pos):.1f}%",
                      ", ".join(gaps[:6]) + ("…" if len(gaps) > 6 else "") or "—"])
    if frows:
        R.table("featcov", "Verification by feature", ["Feature", "Type", "Span", "Verified", "Not verified (1-based)"], frows,
                how="For each annotated feature, the share of its bases verified by at least one read, and the stretches that are not.",
                yours=coverage_yours(rj, comb, focus))
    R.method("Traces", f"Read with Biopython {__import__('Bio').__version__} <code>SeqIO.read(..., 'abi')</code>: channels DATA9–12 "
                       "(FWO_1 dye order), peak locations PLOC2, Phred qualities from <code>letter_annotations['phred_quality']</code>.")
    R.method("Trimming", "Mott's modified trimming algorithm with cutoff 0.05 (the algorithm of Biopython's <code>abi-trim</code> format).")
    R.method("Alignment", "<code>Bio.Align.PairwiseAligner</code>, local mode, match 2, mismatch −3, gap open −6, gap extend −1; "
                          "both strands, the higher-scoring kept" + ("; circular reference aligned as reference + reference" if circ else "") + ".")
    R.method("Calls", "confident: Q ≥ 30 and both neighbours within 2 bases ≥ Q20; check the trace: Q 20–30, a weak neighbour, an "
                      "ambiguous call or a homopolymer context; probably a basecall error: Q < 20. Mixed peaks: second channel ≥ "
                      f"{int(MIXED_RATIO * 100)} % of the called channel. Consequences with the CDS's genetic code (<code>Seq.translate</code>).")
    return R.done()


def diff_label(d):
    if d["type"] == "substitution":
        return f"{d['ref']}{d['ref_pos1']}{d['read']}"
    if d["type"] == "deletion":
        return f"del {d['ref']} at {d['ref_pos1']}"
    return f"ins {d['read']} after {d['ref_pos1'] - 1 if d['ref_pos1'] > 1 else d['ref_pos1']}"


def diff_phrase(d):
    eff = d["effects"][0] if d["effects"] else None
    where = f" in {esc(d['features'][0])}" if d["features"] else ""
    if eff and eff.get("kind") == "missense":
        return f"missense{where} ({eff['hgvs']}, {eff['short']})"
    if eff:
        return f"{esc(eff['text'])}{where}"
    return f"{d['type']} {diff_label(d)}{where}"


def mixed_yours(r):
    m = r["mixed"]
    n = max(1, r["stats"]["trimmed_length"])
    if len(m) > 0.05 * n:
        return (f"{len(m)} positions ({100 * len(m) / n:.1f} % of the trimmed read) have a second peak ≥ {int(MIXED_RATIO * 100)} % — "
                "that is background noise across the read, so single mixed-peak calls here are unreliable.")
    top = sorted(m, key=lambda x: -x["ratio"])[:3]
    return (f"{pl(len(m), 'position')} with a second peak ≥ {int(MIXED_RATIO * 100)} % of the first; the strongest at read "
            + ", ".join(f"{x['read_pos1']} ({x['primary']}/{x['secondary']} {x['ratio']:.2f})" for x in top)
            + ". Look at them in the trace: a clean second peak under a clean first one suggests a mixed colony.")


def _gaps(pos, comb):
    out, start, prev = [], None, None
    for p in pos:
        bad = comb[p] != 2
        if bad and start is None:
            start = p
        if not bad and start is not None:
            out.append(f"{start + 1}" if start == prev else f"{start + 1}–{prev + 1}")
            start = None
        prev = p
    if start is not None:
        out.append(f"{start + 1}" if start == prev else f"{start + 1}–{prev + 1}")
    return out


def coverage_yours(rj, comb, focus):
    L = rj["length"]
    v = sum(1 for x in comb if x == 2)
    diff = sum(1 for x in comb if x == 3)
    s = f"{100 * v / L:.1f} % of {esc(rj['name'])} ({commas(v)} of {commas(L)} bp) is verified at Q ≥ 20"
    if diff:
        s += f"; {commas(diff)} bp carry a difference"
    s += ". "
    targets = [focus] if focus else [f for f in rj["features"] if f["type"] in ("CDS",) or re.search(r"insert", f["name"], re.I)]
    for f in targets[:6]:
        pos = [p for a, b in f["parts"] for p in range(a, b)]
        if not pos:
            continue
        gaps = _gaps(pos, comb)
        fv = sum(1 for p in pos if comb[p] == 2)
        if not gaps:
            s += f"<b>{esc(f['name'])} is fully verified.</b> "
        else:
            s += (f"<b>{esc(f['name'])}</b>: {100 * fv / len(pos):.0f} % verified; not verified: {', '.join(gaps[:4])}"
                  f"{'…' if len(gaps) > 4 else ''}. ")
    if focus is None and targets:
        s += "Pick a feature (the insert) to get a verdict on it alone."
    return s
