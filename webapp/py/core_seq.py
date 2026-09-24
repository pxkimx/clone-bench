"""Construct analysis: reading sequence files, features and CDS checks, restriction analysis, virtual digests
and gels, ORFs, exports, and the small sequence tools. Biopython + numpy + stdlib only."""
from __future__ import annotations

import io
import math
import re
import threading
import warnings

from Bio import BiopythonParserWarning, SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, nt_search
from Bio.SeqUtils.CheckSum import seguid

from .core import Result, UserFacingError, clean_name, commas, esc, file_bytes, pl

DNA_IUPAC = set("ACGTRYSWKMBDHVN")
PROTEIN_ONLY = set("EFILPQ")          # letters that are never nucleotide codes
COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

FORMATS = {  # extension -> (SeqIO format, binary?)
    ".gb": ("genbank", False), ".gbk": ("genbank", False), ".gbff": ("genbank", False), ".genbank": ("genbank", False),
    ".ape": ("genbank", False), ".fa": ("fasta", False), ".fasta": ("fasta", False), ".fna": ("fasta", False),
    ".fas": ("fasta", False), ".fsa": ("fasta", False), ".txt": (None, False), ".seq": (None, False),
    ".embl": ("embl", False), ".dna": ("snapgene", True), ".xdna": ("xdna", True), ".gck": ("gck", True),
}
FORMAT_LABEL = {"genbank": "GenBank", "fasta": "FASTA", "embl": "EMBL", "snapgene": "SnapGene .dna",
                "xdna": "DNA Strider / Serial Cloner .xdna", "gck": "Gene Construction Kit .gck", "raw": "plain sequence"}


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


# ================================================================ reading
def normalize_dna(text: str, what: str = "The sequence") -> tuple[str, list]:
    """Upper-case, strip whitespace/digits, convert RNA, and refuse what is not DNA. Returns (seq, notes)."""
    notes = []
    raw = re.sub(r"[\s\d/]+", "", text or "")
    if not raw:
        raise UserFacingError(f"{what} is empty.")
    lower = sum(1 for c in raw if c.islower())
    s = raw.upper()
    if lower and lower < len(s):
        notes.append(("info", f"{commas(lower)} bases were lower-case (soft-masked or marked). They are analysed like any "
                              "other base; the case is not kept."))
    gaps = sum(s.count(c) for c in "-.~")
    if gaps:
        s = re.sub(r"[-.~]", "", s)
        notes.append(("info", f"{commas(gaps)} gap characters (- or .) were removed."))
    if s.endswith("*"):
        s = s.rstrip("*")
    bad = set(s) - DNA_IUPAC - {"U"}
    if bad:
        punct = [c for c in bad if not c.isalpha()]
        if punct:
            raise UserFacingError(f"{what} does not look like a DNA sequence: it contains ordinary text or punctuation "
                                  f"({' '.join(sorted(punct)[:6])}). Paste only the sequence, or load a GenBank / SnapGene / FASTA file.")
        if bad <= set("EFILPQXZJO"):
            raise UserFacingError(f"{what} looks like a protein sequence (it has {', '.join(sorted(bad)[:6])}). "
                                  "Use the Protein tool for amino-acid sequences.")
        raise UserFacingError(f"{what} contains characters that are not DNA or IUPAC codes: {', '.join(sorted(bad))}.")
    if "U" in s:
        if "T" in s:
            raise UserFacingError(f"{what} contains both U and T. Is it RNA or DNA? Paste it with one or the other.")
        s = s.replace("U", "T")
        notes.append(("warn", "This looked like <b>RNA</b> (U, no T). It was converted to DNA (U → T) for the analysis."))
    n = s.count("N")
    amb = sum(1 for c in s if c not in "ACGTN")
    if n:
        notes.append(("warn", f"{commas(n)} <b>N</b> bases ({100 * n / len(s):.1f}%). Restriction sites, primers and ORFs "
                              "cannot be found inside unknown sequence, so a cutter reported as “single” may cut again in the N stretch."))
    if amb:
        notes.append(("info", f"{commas(amb)} ambiguous IUPAC bases (R, Y, …). They match nothing in site and primer searches."))
    return s, notes


def sniff_format(data: bytes, name: str):
    head = data[:64]
    if head[:4] == b"ABIF":
        raise UserFacingError(f"{name} is a Sanger trace (.ab1), not a construct. Open it in the Sanger tool, "
                              "and load the plasmid map (GenBank / SnapGene / FASTA) here.")
    if data[:1] == b"\t" and data[1:14] == b"\x00\x00\x00\x0eSnapGene":
        return "snapgene", True
    if b"SnapGene" in head[:30]:
        return "snapgene", True
    txt = data[:2000].decode("utf-8", "replace").lstrip()
    if txt.startswith("LOCUS"):
        return "genbank", False
    if txt.startswith(">"):
        return "fasta", False
    if txt.startswith("ID "):
        return "embl", False
    return "raw", False


def parse_sequence_file(f: dict, circular=None) -> tuple[SeqRecord, str, list]:
    """(record, format label, notes) from an uploaded file dict. circular=None keeps the file's own topology."""
    name = f.get("name") or "sequence"
    data = file_bytes(f)
    if not data.strip():
        raise UserFacingError(f"{name} is empty.")
    ext = ("." + name.rsplit(".", 1)[-1].lower()) if "." in name else ""
    if ext in (".ab1", ".abi", ".scf"):
        raise UserFacingError(f"{name} is a Sanger trace, not a construct. Open it in the Sanger tool, and load the "
                              "plasmid map (GenBank / SnapGene / FASTA) here.")
    fmt, binary = FORMATS.get(ext, (None, False))
    if fmt is None:
        fmt, binary = sniff_format(data, name)
    elif data[:4] == b"ABIF":
        sniff_format(data, name)
    notes = []
    if fmt == "raw":
        text = data.decode("utf-8", "replace")
        seq, n2 = normalize_dna(text, f"{name}")
        rec = SeqRecord(Seq(seq), id=re.sub(r"\W+", "_", name.rsplit(".", 1)[0])[:40] or "sequence", name=name, description="")
        rec.annotations["molecule_type"] = "DNA"
        return _finish(rec, circular, "plain sequence", notes + n2, from_text=True)
    handle = io.BytesIO(data) if binary else io.StringIO(data.decode("utf-8", "replace"))
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            recs = list(SeqIO.parse(handle, fmt))
        seen_msgs = set()
        for x in w[:10]:
            msg = str(x.message)
            kind = "wrap" if "origin wrapping" in msg else msg
            if kind in seen_msgs:
                continue
            seen_msgs.add(kind)
            if "origin wrapping" in msg:
                notes.append(("info", "A feature location crossed the origin in a non-standard way (for example "
                                      "<code>3922..2</code>); Biopython read it as a feature that wraps the origin."))
            elif "molecule_type" not in msg:
                notes.append(("info", "Biopython warned while reading the file: " + esc(msg[:200])))
    except Exception as e:  # noqa: BLE001
        hint = ""
        if fmt == "snapgene":
            hint = (" If it was saved by a very old SnapGene (1.x) or a newer format this Biopython does not know, "
                    "export it from SnapGene as GenBank (File → Export → GenBank) and load that.")
        raise UserFacingError(f"{name} could not be read as {FORMAT_LABEL.get(fmt, fmt)}: {str(e)[:200]}.{hint}")
    if not recs:
        raise UserFacingError(f"{name} contains no sequence records. Is it really a {FORMAT_LABEL.get(fmt, fmt)} file?")
    if len(recs) > 1:
        notes.append(("info", f"{name} holds {len(recs)} records; only the first ({esc(recs[0].id)}) was loaded."))
    rec = recs[0]
    try:
        raw = str(rec.seq)
    except UnicodeDecodeError:
        # Biopython keeps sequence as bytes and decodes it as ASCII: a curly quote, Greek letter or other symbol
        # pasted in from Word or a PDF crashed the read instead of being named
        odd = sorted({c for c in bytes(rec.seq).decode("utf-8", "replace") if ord(c) > 127})[:6]
        raise UserFacingError(f"{name}: the sequence contains characters that are not DNA letters ({' '.join(odd)}) — "
                              "usually from copying out of Word or a PDF. Paste it as plain text, or remove them.")
    s, n2 = normalize_dna(raw, name)
    if s != str(rec.seq):
        if len(s) != len(rec) and rec.features:
            raise UserFacingError(f"{name}: the sequence contains gap characters and has annotated features, whose "
                                  "coordinates would no longer fit. Remove the gaps in your editor and load it again.")
        if len(s) != len(rec):
            rec = SeqRecord(Seq(s), id=rec.id, name=rec.name, description=rec.description, annotations=dict(rec.annotations))
        else:
            rec.letter_annotations = {}
            rec.seq = Seq(s)
    return _finish(rec, circular, FORMAT_LABEL.get(fmt, fmt), notes + n2, from_text=fmt == "fasta")


def _finish(rec, circular, label, notes, from_text=False):
    topo = str(rec.annotations.get("topology", "")).lower()
    if circular is None:
        if topo in ("circular", "linear"):
            rec.annotations["topology"] = topo
        else:
            rec.annotations["topology"] = "linear"
            notes.append(("info", f"The {label} file does not say whether the molecule is circular, so it is treated as "
                                  "<b>linear</b>. Tick “circular” and load it again if it is a plasmid."))
    else:
        want = "circular" if circular else "linear"
        if topo and topo != want:
            notes.append(("info", f"The file says <b>{topo}</b>; you asked for <b>{want}</b>, which was used."))
        rec.annotations["topology"] = want
    rec.annotations["molecule_type"] = "DNA"
    return rec, label, notes


def record_from_text(text: str, circular: bool, name: str = "Pasted sequence") -> tuple[SeqRecord, list]:
    text = text or ""
    if text.lstrip().startswith(">"):
        lines = text.strip().splitlines()
        name = lines[0][1:].strip() or name
        body = "".join(l for l in lines[1:] if not l.startswith(">"))
        if sum(1 for l in lines if l.startswith(">")) > 1:
            body = "".join(lines[1:next(i for i, l in enumerate(lines[1:], 1) if l.startswith(">"))])
    else:
        body = text
    s, notes = normalize_dna(body, "The pasted sequence")
    rec = SeqRecord(Seq(s), id=re.sub(r"\W+", "_", name)[:40] or "sequence", name=name[:60], description="")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular" if circular else "linear"
    return rec, notes


# ================================================================ record <-> JSON
def feature_name(f) -> str:
    q = f.qualifiers
    for k in ("label", "gene", "product", "standard_name", "note"):
        if q.get(k):
            v = clean_name(q[k][0])
            return v if k != "note" else v[:40]
    return f.type


def parts_of(loc) -> list:
    return [[int(p.start), int(p.end)] for p in loc.parts]


def feature_span(parts, L, circular):
    """(start0, end0, wraps) of a feature's outer span along the molecule."""
    segs = sorted(parts)
    merged = []
    for s, e in segs:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    if circular and len(merged) > 1 and merged[0][0] == 0 and merged[-1][1] == L:
        # rotate: the arc starts at the last segment and runs through the origin
        return merged[-1][0], merged[0][1], True
    return merged[0][0], merged[-1][1], False


def record_json(rec: SeqRecord) -> dict:
    L = len(rec)
    circ = rec.annotations.get("topology") == "circular"
    feats = []
    for i, f in enumerate(rec.features):
        if f.type == "source" or f.location is None:
            continue
        parts = parts_of(f.location)
        s0, e0, wraps = feature_span(parts, L, circ)
        strand = f.location.strand if f.location.strand in (1, -1) else 0
        q = {k: (v[0] if isinstance(v, list) and v else v) for k, v in f.qualifiers.items()}
        q = {k: str(v)[:4000] for k, v in q.items() if k not in ("ApEinfo_fwdcolor", "ApEinfo_revcolor", "ApEinfo_graphicformat")}
        feats.append({"i": len(feats), "name": feature_name(f), "type": f.type, "strand": strand, "parts": parts,
                      "start0": s0, "end0": e0, "wraps": wraps, "length": len(f.location),
                      "partial5": ("<" in str(f.location.parts[0].start)) if strand >= 0 else (">" in str(f.location.parts[0].end)),
                      "partial3": (">" in str(f.location.parts[-1].end)) if strand >= 0 else ("<" in str(f.location.parts[-1].start)),
                      "fuzzy": any(("<" in str(p.start)) or (">" in str(p.end)) for p in f.location.parts),
                      "qualifiers": q})
    return {"name": rec.name if rec.name and rec.name != "<unknown name>" else rec.id, "id": rec.id,
            "description": rec.description if rec.description != "<unknown description>" else "",
            "length": L, "circular": circ, "seq": str(rec.seq), "features": feats}


def record_from_json(rj: dict) -> SeqRecord:
    if not isinstance(rj, dict) or not rj.get("seq"):
        raise UserFacingError("No construct was given. Load one in the Construct tool first.")
    s, _ = normalize_dna(rj["seq"], "The construct")
    rec = SeqRecord(Seq(s), id=re.sub(r"\W+", "_", rj.get("id") or rj.get("name") or "construct")[:40],
                    name=(rj.get("name") or "construct")[:60], description=rj.get("description") or "")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular" if rj.get("circular") else "linear"
    for fj in rj.get("features") or []:
        loc = feature_location(fj)
        if loc is None:
            continue
        q = {k: [v] for k, v in (fj.get("qualifiers") or {}).items()}
        q.setdefault("label", [fj.get("name") or fj.get("type")])
        rec.features.append(SeqFeature(loc, type=fj.get("type") or "misc_feature", qualifiers=q))
    return rec


def feature_location(fj):
    from Bio.SeqFeature import AfterPosition, BeforePosition
    strand = fj.get("strand") or None
    locs = [SimpleLocation(int(s), int(e), strand=strand) for s, e in fj.get("parts") or [] if int(e) > int(s)]
    if locs and fj.get("partial3"):       # e.g. an ORF that runs off the end of the sequence: 34..>1074
        last = locs[-1]
        locs[-1] = SimpleLocation(last.start, AfterPosition(int(last.end)), strand=strand) if strand != -1 else \
            SimpleLocation(BeforePosition(int(last.start)), last.end, strand=strand)
    if locs and fj.get("partial5"):
        first = locs[0]
        locs[0] = SimpleLocation(BeforePosition(int(first.start)), first.end, strand=strand) if strand != -1 else \
            SimpleLocation(first.start, AfterPosition(int(first.end)), strand=strand)
    if not locs:
        return None
    return locs[0] if len(locs) == 1 else CompoundLocation(locs)


def extract_feature(seq: str, fj: dict) -> str:
    out = "".join(seq[s:e] for s, e in fj["parts"]) if fj.get("strand", 1) != -1 else \
        "".join(revcomp(seq[s:e]) for s, e in fj["parts"])
    return out


def find_feature(rj: dict, key):
    """A feature by index, or by name (case-insensitive, exact first, then substring)."""
    feats = rj.get("features") or []
    if key is None or key == "":
        return None
    if isinstance(key, int) or (isinstance(key, str) and key.isdigit()):
        k = int(key)
        return feats[k] if 0 <= k < len(feats) else None
    k = str(key).lower()
    for f in feats:
        if f["name"].lower() == k:
            return f
    for f in feats:
        if k in f["name"].lower():
            return f
    return None


def fmt_span(fj, L=None):
    s, e = fj["start0"] + 1, fj["end0"]
    return f"{commas(s)}–{commas(e)}" + (" (across origin)" if fj.get("wraps") else "")


# ================================================================ CDS checks
def check_cds(fj: dict, seq: str) -> dict:
    """Translate a CDS feature and compare it with its /translation. Returns {translation, issues, ok, notes}."""
    q = fj.get("qualifiers") or {}
    nt = extract_feature(seq, fj)
    try:
        table = int(q.get("transl_table") or 11)
        ct = CodonTable.unambiguous_dna_by_id[table]
    except Exception:  # noqa: BLE001
        table, ct = 11, CodonTable.unambiguous_dna_by_id[11]
    try:
        cs = int(q.get("codon_start") or 1)
    except ValueError:
        cs = 1
    nt = nt[cs - 1:]
    issues, notes = [], []
    if len(nt) % 3:
        issues.append(f"length {commas(len(nt))} bp is not a multiple of 3")
    first = nt[:3]
    starts_ok = first in ct.start_codons
    if not starts_ok and not fj.get("fuzzy"):
        issues.append(f"does not begin with a start codon (begins {first})")
    elif starts_ok and first != "ATG":
        notes.append(f"starts with {first}, a valid alternative start in NCBI table {table}")
    usable = nt[:len(nt) // 3 * 3]
    aa = str(Seq(usable).translate(table=table)) if usable else ""
    if aa and starts_ok:
        aa = "M" + aa[1:]
    has_stop = aa.endswith("*")
    if not has_stop and not fj.get("fuzzy"):
        issues.append("does not end with a stop codon")
    body = aa[:-1] if has_stop else aa
    internal = [i + 1 for i, c in enumerate(body) if c == "*"]
    if internal and not has_stop and len(internal) == 1 and internal[0] >= len(body) - 3:
        # the usual annotation slip: the feature runs a codon or two past the real stop
        issues[-1] = (f"the first in-frame stop is codon {internal[0]}, {3 * (len(body) - internal[0])} bp before the annotated "
                      "end — the feature runs past its stop codon")
    elif internal:
        issues.append(f"internal stop codon{'s' if len(internal) > 1 else ''} at residue {', '.join(map(str, internal[:5]))}")
    calc = body.split("*")[0] if internal and not has_stop and len(internal) == 1 and internal[0] >= len(body) - 3 else body
    given = q.get("translation")
    match = None
    if given:
        g = re.sub(r"[\s,]", "", given).rstrip("*")
        if "," in given:
            notes.append("the file's /translation contains commas where SnapGene joins feature segments; they were ignored")
        if g == calc:
            match = True
        else:
            match = False
            k = next((i for i, (a, b) in enumerate(zip(g, calc)) if a != b), min(len(g), len(calc)))
            detail = (f"first difference at residue {k + 1}: file has {g[k] if k < len(g) else '(end)'}, "
                      f"the sequence translates to {calc[k] if k < len(calc) else '(end)'}")
            if len(g) != len(calc):
                detail += f"; lengths {len(g)} vs {len(calc)} aa"
            issues.append("the translation in the file does not match the DNA (" + detail + ")")
    if "X" in calc:
        notes.append(f"{calc.count('X')} codons contain ambiguous bases and translate to X")
    return {"translation": calc, "table": table, "issues": issues, "notes": notes, "match": match, "aa": len(calc)}


# ================================================================ restriction analysis
_RLOCK = threading.Lock()   # Bio.Restriction keeps search state on the enzyme classes: one search at a time
_BATCH = None
_INFO = {}


# Commercially available, but they cut only DNA carrying 5-methyl- or 5-hydroxymethylcytosine (AbaSI needs
# glucosylated 5hmC). A plasmid prep has none, so they never cut it — yet searched as plain sequence patterns
# they "cut" at nearly every C (AbaSI's site is one base) and filled the enzyme table with nonsense.
MODIFICATION_DEPENDENT = {"AbaSI", "FspEI", "LpnPI", "MspJI", "SgeI"}


def _batch():
    global _BATCH
    from Bio.Restriction import CommOnly, RestrictionBatch
    if _BATCH is None:
        _BATCH = RestrictionBatch([e for e in CommOnly if str(e) not in MODIFICATION_DEPENDENT])
    return _BATCH


def enzyme(name: str):
    from Bio.Restriction import AllEnzymes
    for e in AllEnzymes:
        if str(e).lower() == str(name).strip().lower():
            return e
    raise UserFacingError(f"“{name}” is not a restriction enzyme Biopython knows. Check the spelling (e.g. EcoRI, BamHI, XhoI).")


def enzyme_info(e) -> dict:
    k = str(e)
    if k in _INFO:
        return _INFO[k]
    site = str(e.site)
    if e.is_blunt():
        oh = "blunt"
    elif e.is_5overhang():
        oh = f"5′ {abs(e.ovhg)} nt"
    elif e.is_3overhang():
        oh = f"3′ {abs(e.ovhg)} nt"
    else:
        oh = "unknown"
    info = {"name": k, "site": site, "elucidate": e.elucidate(), "overhang": oh, "size": e.size,
            "degenerate": bool(set(site) - set("ACGT")), "outside": bool(e.fst5 is not None and (e.fst5 < 0 or e.fst5 > e.size)),
            "suppliers": "".join(sorted(e.suppl)), "neb": "N" in e.suppl}
    _INFO[k] = info
    return info


def cut_map(seq: str, circular: bool, names=None) -> dict:
    """{enzyme name: sorted cut positions} — Biopython's convention: the 1-based position of the first base
    after the cut on the top strand. names=None searches every commercially available enzyme (CommOnly)."""
    from Bio.Restriction import Analysis, RestrictionBatch
    with _RLOCK:
        batch = _batch() if names is None else RestrictionBatch([enzyme(n) for n in names])
        an = Analysis(batch, Seq(seq), linear=not circular)
        m = {str(e): sorted(set(int(x) for x in v)) for e, v in an.mapping.items()}
    return m


def features_at(rj: dict, pos0: int, types_skip=("source",)) -> list:
    """Names of features that contain a cut made just before 0-based base pos0 (strictly inside, not at an edge)."""
    out = []
    for f in rj["features"]:
        if f["type"] in types_skip:
            continue
        for s, e in f["parts"]:
            if s < pos0 < e:
                out.append(f["name"])
                break
    return out


def in_span(fj, p0, L):
    """Is 0-based cut coordinate p0 (cut between p0-1 and p0) strictly inside the feature's outer span?"""
    s, e, w = fj["start0"], fj["end0"], fj["wraps"]
    if not w:
        return s < p0 < e
    return p0 > s or p0 < e


def enzyme_rows(rj: dict, mapping: dict, only=None) -> list:
    L = rj["length"]
    rows = []
    for name, cuts in mapping.items():
        if only is not None and name not in only:
            continue
        if not cuts:
            continue
        info = enzyme_info(enzyme(name))
        feats = []
        for p in cuts:
            for n in features_at(rj, (p - 1) % L):
                if n not in feats:
                    feats.append(n)
        rows.append({"enzyme": name, "cuts": len(cuts), "positions": cuts, "site": info["elucidate"],
                     "overhang": info["overhang"], "features": feats, "neb": info["neb"], "size": info["size"],
                     "degenerate": info["degenerate"], "suppliers": info["suppliers"]})
    rows.sort(key=lambda r: (r["cuts"], r["enzyme"].lower()))
    return rows


def filter_enzymes(rj: dict, mapping: dict, flt: str, feature=None) -> dict:
    """Named filters. 'outside' and 'inside' use Bio.Restriction.Analysis.do_not_cut / only_between on the
    selected feature (by hand when the feature wraps the origin, which Analysis's boundaries cannot express)."""
    L, circ = rj["length"], rj["circular"]
    fj = find_feature(rj, feature) if feature not in (None, "") else None
    if flt in ("outside", "inside", "flank") and fj is None:
        raise UserFacingError("Pick a feature first (click it on the map or choose it from the list).")
    single = {k: v for k, v in mapping.items() if len(v) == 1}
    note = ""
    if flt == "single":
        keep = single
    elif flt == "double":
        keep = {k: v for k, v in mapping.items() if len(v) == 2}
    elif flt == "none":
        keep = {k: v for k, v in mapping.items() if not v}
    elif flt == "all":
        keep = {k: v for k, v in mapping.items() if v}
    elif flt in ("outside", "flank", "inside"):
        if not fj["wraps"]:
            from Bio.Restriction import Analysis, RestrictionBatch
            names = list(mapping)
            with _RLOCK:
                an = Analysis(RestrictionBatch([enzyme(n) for n in names]), Seq(rj["seq"]), linear=not circ)
                # a cut strictly inside [start0, end0) has position p (first base after the cut, 1-based) with
                # start0 + 2 <= p <= end0, i.e. Analysis's half-open test start <= p < end with these borders
                lo, hi = fj["start0"] + 2, fj["end0"] + 1
                if flt == "inside":
                    d = an.only_between(lo, hi) if hi > lo else {}
                else:
                    d = an.do_not_cut(lo, hi) if hi > lo else dict(an.mapping)
                d = {str(k): sorted(v) for k, v in d.items()}
            if flt == "inside":
                keep = {k: v for k, v in d.items() if v}
            else:
                keep = {k: v for k, v in d.items() if len(v) == 1}
            note = f"Computed with Analysis.{'only_between' if flt == 'inside' else 'do_not_cut'}({lo}, {hi})."
        else:
            inside = lambda p: in_span(fj, (p - 1) % L, L)  # noqa: E731
            if flt == "inside":
                keep = {k: v for k, v in mapping.items() if v and all(inside(p) for p in v)}
            else:
                keep = {k: v for k, v in single.items() if not inside(v[0])}
            note = "The feature wraps the origin, so the inside/outside test was done by hand on both arcs."
    else:
        raise UserFacingError(f"Unknown filter '{flt}'.")
    out = {"filter": flt, "feature": fj["name"] if fj else None, "note": note}
    if flt == "none":
        out["names"] = sorted(keep, key=str.lower)
        out["rows"] = []
        return out
    rows = enzyme_rows(rj, keep)
    if flt == "flank" and fj:
        s, e = fj["start0"], fj["end0"]
        up, down = [], []
        for r in rows:
            p0 = (r["positions"][0] - 1) % L
            du = (s - p0) % L if circ else s - p0
            dd = (p0 - e) % L if circ else p0 - e
            if not circ and du < 0 and dd < 0:
                continue
            if circ:
                (up if du <= dd else down).append((du if du <= dd else dd, r))
            else:
                (up if du >= 0 else down).append((du if du >= 0 else dd, r))
        up.sort(key=lambda x: x[0])
        down.sort(key=lambda x: x[0])
        out["upstream"] = [dict(r, distance=d) for d, r in up[:15]]
        out["downstream"] = [dict(r, distance=d) for d, r in down[:15]]
    out["rows"] = rows
    return out


def single_cutter_labels(rj: dict, mapping: dict) -> list:
    """Single cutters grouped by cut position (isoschizomers share a label)."""
    groups = {}
    for name, cuts in mapping.items():
        if len(cuts) != 1:
            continue
        info = enzyme_info(enzyme(name))
        groups.setdefault(cuts[0], []).append(info)
    out = []
    for pos, infos in sorted(groups.items()):
        infos.sort(key=lambda i: (not i["neb"], i["degenerate"], -i["size"], len(i["name"]), i["name"]))
        out.append({"pos": pos, "names": [i["name"] for i in infos], "neb": any(i["neb"] for i in infos),
                    "six": any(i["size"] >= 6 and not i["degenerate"] for i in infos),
                    "sites": sorted({i["elucidate"] for i in infos})})
    return out


# ================================================================ digest & gel
LADDERS = {
    "generuler_1kb_plus": {"label": "GeneRuler 1 kb Plus (Thermo SM1331)",
                           "bands": [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75]},
    "neb_1kb_plus": {"label": "NEB 1 kb Plus (N3200)",
                     "bands": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500,
                               400, 300, 200, 100]},
}
# enzymes most labs have in the freezer: preferred for default lanes and suggested digests
COMMON = ["EcoRI", "BamHI", "HindIII", "XhoI", "NotI", "XbaI", "KpnI", "SacI", "PstI", "SalI", "NcoI", "NdeI", "SpeI",
          "BglII", "EcoRV", "NheI", "ApaI", "SmaI", "ClaI", "SphI", "AgeI", "MluI", "NsiI", "PvuI", "PvuII", "ScaI", "AflII",
          "AvrII", "BsrGI", "PacI", "AscI", "SbfI", "NruI", "StuI", "HpaI", "SacII", "MfeI", "BspHI", "AatII", "SnaBI"]
AGAROSE = {0.7: (800, 12000), 1.0: (500, 10000), 2.0: (100, 2000)}   # well-resolved size range (bp)


def migration(size: float, pct: float) -> tuple[float, str]:
    """Approximate band position (0 = well, 1 = gel bottom) and a note. Log-linear inside the resolvable range."""
    lo, hi = AGAROSE.get(float(pct), AGAROSE[1.0])
    a, b = 0.08, 0.90
    if size >= hi:
        y = a * (hi / size) ** 0.35       # compresses toward the well
        return max(0.025, y), "above the resolvable range: compressed near the well"
    if size < lo:
        y = b + (1 - b) * min(1.0, math.log(lo / size) / math.log(4))
        if size < lo / 2:
            return min(0.99, y), "below the resolvable range: diffuse, may run off"
        return y, "below the well-resolved range"
    return a + (b - a) * (math.log(hi) - math.log(size)) / (math.log(hi) - math.log(lo)), ""


def fragments(L: int, circular: bool, cuts: dict) -> list:
    """cuts: {pos0: [enzyme names]} with pos0 = 0-based index of the first base after the cut."""
    pos = sorted(cuts)
    if not pos:
        return [{"size": L, "start0": 0, "end0": L, "left": "end" if not circular else None,
                 "right": "end" if not circular else None, "uncut": True}]
    out = []
    if circular:
        for i, p in enumerate(pos):
            q = pos[(i + 1) % len(pos)]
            size = (q - p) % L or L
            out.append({"size": size, "start0": p, "end0": q, "left": " + ".join(cuts[p]), "right": " + ".join(cuts[q]),
                        "wraps": q <= p})
    else:
        bounds = [0] + [p for p in pos if 0 < p < L] + [L]
        for a, b in zip(bounds, bounds[1:]):
            if b > a:
                out.append({"size": b - a, "start0": a, "end0": b, "left": " + ".join(cuts[a]) if a in cuts else "end",
                            "right": " + ".join(cuts[b]) if b in cuts else "end"})
    return out


def bands_for(frags: list, pct: float) -> list:
    """Group fragments that co-migrate (within ~5 % of each other's size) into bands."""
    fs = sorted(frags, key=lambda f: -f["size"])
    groups = []
    for f in fs:
        if groups and (groups[-1][-1]["size"] - f["size"]) / groups[-1][-1]["size"] <= 0.05:
            groups[-1].append(f)
        else:
            groups.append([f])
    out = []
    for g in groups:
        size = sum(f["size"] for f in g) / len(g)
        y, note = migration(size, pct)
        mass = sum(f["size"] for f in g)
        out.append({"size": round(size), "sizes": [f["size"] for f in g], "n": len(g), "y": y, "mass": mass, "note": note,
                    "faint": size < 100})
    return out


def digest(rj: dict, enzymes: list, pct: float = 1.0) -> dict:
    L, circ = rj["length"], rj["circular"]
    enzymes = [e for e in (enzymes or []) if str(e).strip()]
    if len(enzymes) > 3:
        raise UserFacingError("Pick at most three enzymes per lane.")
    if not enzymes:
        uncut = True
        frs = fragments(L, circ, {})
    else:
        m = cut_map(rj["seq"], circ, [str(e).strip() for e in enzymes])
        cuts = {}
        for n, ps in m.items():
            for p in ps:
                cuts.setdefault((p - 1) % L, []).append(n)
        frs = fragments(L, circ, cuts)
        uncut = not cuts
    if uncut and circ:
        return {"enzymes": enzymes, "uncut_circular": True, "fragments": [], "cuts": 0,
                "bands": [{"size": None, "y": migration(L * 0.65, pct)[0], "mass": L, "n": 1, "supercoiled": True,
                           "note": "supercoiled: runs faster than its size", "sizes": []}]}
    for f in frs:
        f["start"] = f["start0"] + 1
        f["end"] = f["end0"] if f["end0"] else L
    return {"enzymes": enzymes, "fragments": frs, "cuts": 0 if uncut else len(frs) - (0 if circ else 1),
            "bands": bands_for(frs, pct)}


def without_feature(rj: dict, fj: dict) -> dict:
    """The construct with one feature's span deleted — a stand-in for 'the empty vector'."""
    seq = rj["seq"]
    s, e = fj["start0"], fj["end0"]
    if fj["wraps"]:
        new = seq[e:s]
    else:
        new = seq[:s] + seq[e:]
    return {"name": rj["name"] + " minus " + fj["name"], "length": len(new), "circular": rj["circular"], "seq": new,
            "features": []}


def lane_summary(lane, pct):
    fr = lane.get("fragments") or []
    if lane.get("uncut_circular"):
        return "uncut circular DNA — supercoiled, nicked and linear forms run at different, size-independent positions"
    bands = lane["bands"]
    co = [b for b in bands if b["n"] > 1]
    faint = [f for f in fr if f["size"] < 100]
    lo, hi = AGAROSE.get(float(pct), AGAROSE[1.0])
    out_hi = [f for f in fr if f["size"] >= hi]
    parts = [f"{pl(len(fr), 'fragment')} → {pl(len(bands), 'visible band')}"]
    if co:
        parts.append("; ".join(f"{' + '.join(commas(s) for s in b['sizes'])} bp co-migrate as one brighter band" for b in co))
    if faint:
        parts.append(f"{' and '.join(commas(f['size']) for f in faint)} bp {'is' if len(faint) == 1 else 'are'} under 100 bp: faint or off the gel")
    if out_hi:
        parts.append(f"{' and '.join(commas(f['size']) for f in out_hi)} bp {'runs' if len(out_hi) == 1 else 'run'} above the range a {pct:g} % gel resolves")
    return ". ".join(parts)


def band_signature(bands):
    return [b["size"] for b in bands if b.get("size") and b["size"] >= 100]


def patterns_differ(b1, b2) -> bool:
    s1, s2 = band_signature(b1), band_signature(b2)
    if len(s1) != len(s2):
        return True
    for x in s1:
        if not any(abs(x - y) / max(x, y) <= 0.05 for y in s2):
            return True
    return False


def digest_api(args: dict) -> dict:
    """Lanes of a virtual gel. args: record, lanes [[enzymes]], ladder, agarose, compare_feature (optional)."""
    rj = args.get("record")
    if not rj:
        raise UserFacingError("Load a construct first.")
    pct = float(args.get("agarose") or 1.0)
    if pct not in AGAROSE:
        pct = min(AGAROSE, key=lambda k: abs(k - pct))
    ladder_key = args.get("ladder") or "generuler_1kb_plus"
    ladder = LADDERS.get(ladder_key, LADDERS["generuler_1kb_plus"])
    lanes = []
    for enz in (args.get("lanes") or [[]])[:8]:
        d = digest(rj, enz, pct)
        d["label"] = " + ".join(enz) if enz else "uncut"
        d["summary"] = lane_summary(d, pct)
        lanes.append(d)
    cmp_lanes = []
    fj = find_feature(rj, args.get("compare_feature")) if args.get("compare_feature") not in (None, "") else None
    if fj:
        emp = without_feature(rj, fj)
        for enz, d0 in zip(args.get("lanes") or [], lanes):
            d = digest(emp, enz, pct)
            d["label"] = (" + ".join(enz) if enz else "uncut") + " · no " + fj["name"]
            d["summary"] = lane_summary(d, pct)
            d["empty"] = True
            d["differs"] = patterns_differ(d0["bands"], d["bands"]) if not d0.get("uncut_circular") else False
            cmp_lanes.append(d)
    ladder_bands = []
    for s in ladder["bands"]:
        y, note = migration(s, pct)
        ladder_bands.append({"size": s, "y": y, "note": note})
    # the "In your data" text for the gel
    yours = []
    for d in lanes:
        yours.append(f"<b>{esc(d['label'])}</b>: {d['summary']}.")
    for d in cmp_lanes:
        if d.get("uncut_circular"):
            yours.append(f"<b>{esc(d['label'])}</b>: uncut plasmids are not compared — supercoiled DNA does not run at its size.")
            continue
        yours.append(f"<b>{esc(d['label'])}</b>: {d['summary']}. "
                     + ("The pattern <b>differs</b> from your construct, so this digest tells them apart."
                        if d["differs"] else "The pattern is <b>the same within gel resolution</b> — this digest does not tell them apart."))
    verdict = []
    for d in lanes:
        if d.get("uncut_circular") or not d["fragments"]:
            continue
        clean = all(b["n"] == 1 for b in d["bands"]) and all(100 <= f["size"] < AGAROSE[pct][1] for f in d["fragments"])
        if len(d["fragments"]) == 1:
            verdict.append(f"{esc(d['label'])} only linearises the plasmid: it confirms the total size "
                           f"({commas(d['fragments'][0]['size'])} bp), not the arrangement inside it.")
        elif clean:
            verdict.append(f"{esc(d['label'])} gives {len(d['fragments'])} separate, well-resolved bands on a {pct:g} % gel — a clean fingerprint of this construct.")
        else:
            verdict.append(f"{esc(d['label'])} does not give every fragment its own resolvable band on a {pct:g} % gel.")
    lo, hi = AGAROSE[pct]
    return {"lanes": lanes, "compare": cmp_lanes, "ladder": {"key": ladder_key, "label": ladder["label"], "bands": ladder_bands},
            "agarose": pct, "range": [lo, hi], "yours": " ".join(yours + verdict),
            "compare_feature": fj["name"] if fj else None}


def suggest_digests(rj: dict, feature=None, pct: float = 1.0, top: int = 5) -> list:
    """Enzyme combinations (one or two NEB-sold enzymes with 6+ bp unambiguous sites) whose fragments are all
    resolvable, well spread, and — when a feature is given — cut inside it and differ from the construct
    without it. A search over sensible candidates, scored by how far apart the closest two bands sit."""
    L, circ = rj["length"], rj["circular"]
    m = cut_map(rj["seq"], circ)
    fj = find_feature(rj, feature) if feature not in (None, "") else None
    cand = {}
    for n, ps in m.items():
        if not 1 <= len(ps) <= 3:
            continue
        info = enzyme_info(enzyme(n))
        if info["size"] < 6 or info["degenerate"] or not info["neb"] or info["outside"]:
            continue
        key = tuple(ps)
        if key not in cand or len(n) < len(cand[key]):
            cand[key] = n      # isoschizomers cut identically: keep one
    names = list(cand.values())
    lo, hi = AGAROSE[pct]
    lo = max(lo, 250)
    emp = without_feature(rj, fj) if fj else None
    em = cut_map(emp["seq"], circ, names) if emp else None
    combos = [[n] for n in names if len(m[n]) >= 2] + \
             [[a, b] for i, a in enumerate(names) for b in names[i + 1:] if len(m[a]) + len(m[b]) <= 4]
    scored = []
    for combo in combos:
        cuts = {}
        for n in combo:
            for p in m[n]:
                cuts.setdefault((p - 1) % L, []).append(n)
        if len(cuts) < 2:
            continue
        frs = fragments(L, circ, cuts)
        sizes = sorted(f["size"] for f in frs)
        if sizes[0] < lo or sizes[-1] > hi or len(sizes) > 4:
            continue
        gaps = [math.log(b / a) for a, b in zip(sizes, sizes[1:])]
        spread = min(gaps) if gaps else 0
        if spread < math.log(1.10):
            continue
        why = []
        score = min(spread, math.log(2.0)) + 0.25 * (len(sizes) - 2) + 0.15 * sum(1 for n in combo if n in COMMON)
        if fj:
            inside = [p for p in cuts if in_span(fj, p, L)]
            if not inside:
                continue
            ecuts = {}
            for n in combo:
                for p in em[n]:
                    ecuts.setdefault((p - 1) % len(emp["seq"]), []).append(n)
            ef = fragments(len(emp["seq"]), circ, ecuts)
            if not patterns_differ(bands_for(frs, pct), bands_for(ef, pct)):
                continue
            why.append(f"cuts inside {fj['name']}; without it you would see {' + '.join(commas(f['size']) for f in sorted(ef, key=lambda f: -f['size']))} bp")
            score += 0.5
        if len(combo) == 1:
            score += 0.05     # one enzyme is simpler than two
        scored.append((score, combo, sorted(sizes, reverse=True), why))
    scored.sort(key=lambda x: -x[0])
    return [{"enzymes": c, "sizes": s, "min_ratio": round(min(b / a for a, b in zip(sorted(s), sorted(s)[1:])), 2),
             "why": "; ".join(w)} for sc, c, s, w in scored[:top]]


def suggest_api(args):
    rj = args.get("record")
    if not rj:
        raise UserFacingError("Load a construct first.")
    pct = float(args.get("agarose") or 1.0)
    s = suggest_digests(rj, args.get("feature"), pct if pct in AGAROSE else 1.0)
    return {"suggestions": s, "feature": args.get("feature"),
            "note": "Searched one or two NEB-sold enzymes with 6+ bp unambiguous sites; every fragment between "
                    f"250 bp and {commas(AGAROSE.get(pct, AGAROSE[1.0])[1])} bp and neighbouring bands at least 10 % apart."}


def enzymes_api(args):
    rj = args.get("record")
    if not rj:
        raise UserFacingError("Load a construct first.")
    m = cut_map(rj["seq"], rj["circular"])
    return filter_enzymes(rj, m, args.get("filter") or "single", args.get("feature"))


# ================================================================ ORFs
def codon_tables_api(args=None):
    return {"tables": [{"id": i, "name": t.names[0] if t.names else str(i)}
                       for i, t in sorted(CodonTable.unambiguous_dna_by_id.items())]}


def find_orfs(seq: str, circular: bool, table: int = 11, min_aa: int = 75, starts: str = "atg") -> list:
    """Six-frame ORF scan: from the first start codon after an in-frame stop to the next stop.

    A circular sequence is scanned as sequence + sequence, so an ORF can cross the origin; each ORF is kept once
    (the longest one ending at a given stop codon), in the coordinates of the original sequence."""
    L = len(seq)
    try:
        ct = CodonTable.unambiguous_dna_by_id[int(table)]
    except (KeyError, ValueError):
        raise UserFacingError(f"NCBI genetic code {table} does not exist.")
    stops = set(ct.stop_codons)
    start_set = {"ATG"} if starts == "atg" else set(ct.start_codons)
    best = {}
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        ext = s + s if circular else s
        n = len(ext)
        for frame in range(3):
            open_at = None
            for i in range(frame, n - 2, 3):
                if i >= L and open_at is None:
                    break                      # anything starting here is a copy of an ORF already seen
                cod = ext[i:i + 3]
                if cod in stops:
                    if open_at is not None:
                        a, b = open_at, i + 3
                        if b - a <= L and (b - a - 3) // 3 >= min_aa:
                            key = (strand, b % L if circular else b)
                            if key not in best or best[key][1] - best[key][0] < b - a:
                                best[key] = (a, b)
                        open_at = None
                elif open_at is None and cod in start_set:
                    open_at = i
    res = []
    for (strand, _), (a, b) in best.items():
        s = seq if strand == 1 else revcomp(seq)
        ext = s + s if circular else s
        prot = str(Seq(ext[a:b]).translate(table=int(table)))
        prot = ("M" + prot[1:]).rstrip("*")
        wraps = b > L
        if strand == 1:
            s0, e0 = a, (b - 1) % L + 1
            parts = [[s0, L], [0, e0]] if wraps else [[s0, e0]]
        else:
            s0, e0 = (L - b) % L, L - a
            parts = [[0, e0], [s0, L]] if wraps else [[s0, e0]]
        frame = a % 3 if not circular else a % 3
        res.append({"strand": strand, "frame": ("+" if strand == 1 else "−") + str(frame + 1), "start0": s0, "end0": e0,
                    "wraps": wraps, "parts": parts, "aa": len(prot), "nt": b - a, "protein": prot})
    res.sort(key=lambda r: (-r["aa"], r["start0"]))
    return res


def match_orfs(rj, orfs):
    cds = [f for f in rj["features"] if f["type"] == "CDS"]
    L = rj["length"]
    for o in orfs:
        o["match"] = ""
        stop_end = o["end0"] if o["strand"] == 1 else o["start0"]
        for f in cds:
            if f["strand"] != o["strand"]:
                continue
            fstop = f["end0"] if f["strand"] == 1 else f["start0"]
            fstart = f["start0"] if f["strand"] == 1 else f["end0"]
            ostart = o["start0"] if o["strand"] == 1 else o["end0"]
            if fstop % L == stop_end % L:
                o["match"] = f"= {f['name']}" if fstart % L == ostart % L else f"same stop as {f['name']} (different start)"
                break
        if not o["match"]:
            for f in cds:
                fstart = f["start0"] if f["strand"] == 1 else f["end0"]
                ostart = o["start0"] if o["strand"] == 1 else o["end0"]
                if f["strand"] == o["strand"] and fstart % L == ostart % L:
                    o["match"] = f"same start as {f['name']} (different stop)"
                    break
        if not o["match"]:
            for f in cds:
                if spans_overlap(o, f, L):
                    o["match"] = f"overlaps {f['name']}"
                    break
    return orfs


def spans_overlap(a, b, L):
    def segs(x):
        return x["parts"] if "parts" in x else [[x["start0"], x["end0"]]]
    for s1, e1 in segs(a):
        for s2, e2 in segs(b):
            if s1 < e2 and s2 < e1:
                return True
    return False


def orfs_api(args):
    rj = args.get("record")
    if not rj:
        raise UserFacingError("Load a construct first.")
    table = int(args.get("table") or 11)
    min_aa = max(10, int(args.get("min_aa") or 75))
    starts = args.get("starts") or "atg"
    orfs = match_orfs(rj, find_orfs(rj["seq"], rj["circular"], table, min_aa, starts))
    return {"orfs": orfs[:300], "n": len(orfs), "table": table, "min_aa": min_aa, "starts": starts,
            "yours": orf_yours(rj, orfs, table, min_aa, starts)}


def orf_yours(rj, orfs, table, min_aa, starts):
    if not orfs:
        return f"No open reading frame of {min_aa} codons or more in any of the six frames (NCBI table {table})."
    matched = [o for o in orfs if o["match"].startswith("=")]
    novel = [o for o in orfs if not o["match"]]
    cds = [f for f in rj["features"] if f["type"] == "CDS"]
    wrap = [o for o in orfs if o["wraps"]]
    o0 = orfs[0]
    t = (f"{pl(len(orfs), 'ORF')} of ≥ {min_aa} codons ({'ATG' if starts == 'atg' else 'any table-' + str(table) + ' start'} to stop). "
         f"The longest is {o0['aa']} aa on the {'forward' if o0['strand'] == 1 else 'reverse'} strand, "
         f"{commas(o0['start0'] + 1)}–{commas(o0['end0'])}{' across the origin' if o0['wraps'] else ''}"
         f"{' (' + o0['match'] + ')' if o0['match'] else ''}. ")
    if cds:
        t += f"{len(matched)} of {len(cds)} annotated CDS are recovered exactly. "
    if novel:
        t += (f"{pl(len(novel), 'ORF')} overlap{'s' if len(novel) == 1 else ''} no annotated CDS — long ORFs on the antisense strand of a real gene are common and "
              "usually not expressed, so treat these as candidates, not genes. ")
    if wrap:
        t += f"{pl(len(wrap), 'ORF')} cross{'es' if len(wrap) == 1 else ''} the origin."
    return t


# ================================================================ construct (the whole tool)
def construct(args: dict) -> dict:
    circ_arg = args.get("circular")
    circ_arg = None if circ_arg in (None, "", "auto") else bool(circ_arg)
    if args.get("file"):
        rec, label, notes = parse_sequence_file(args["file"], circ_arg)
        name = args["file"].get("name") or rec.id
    elif args.get("text"):
        rec, notes = record_from_text(args["text"], bool(circ_arg), args.get("name") or "Pasted sequence")
        label, name = "pasted sequence", args.get("name") or rec.name
    else:
        raise UserFacingError("Drop a sequence file or paste a sequence.")
    rj = record_json(rec)
    L, circ, seq = rj["length"], rj["circular"], rj["seq"]
    if L < 20:
        raise UserFacingError(f"The sequence is only {L} bp — too short to be a construct. For a primer, use the Primers tool.")
    R = Result("construct", name)
    R.r["record"] = rj
    gc = 100 * gc_fraction(seq, ambiguous="ignore")

    # ---- CDS checks
    cds_rows, n_bad, n_ok = [], 0, 0
    for f in rj["features"]:
        if f["type"] != "CDS":
            continue
        chk = check_cds(f, seq)
        f["translation"] = chk["translation"]
        f["check"] = chk
        if chk["issues"]:
            n_bad += 1
            R.flag("warn", f"CDS <b>{esc(f['name'])}</b> ({fmt_span(f)}): " + "; ".join(chk["issues"]) + ".")
        else:
            n_ok += 1
        cds_rows.append(f)
    if n_ok:
        checked = [f for f in cds_rows if not f["check"]["issues"]]
        with_tr = [f for f in checked if f["check"]["match"]]
        R.flag("ok", f"{pl(n_ok, 'CDS', 'CDS')} translate cleanly (start codon, one stop at the end, length a multiple of 3)"
                     + (f"; {len(with_tr)} of them match the <code>/translation</code> stored in the file" if with_tr else "")
                     + ": " + ", ".join(esc(f["name"]) for f in checked[:8]) + ".")
    for lvl, t in notes:
        R.flag(lvl, t)
    wraps = [f for f in rj["features"] if f["wraps"]]
    if wraps:
        R.flag("info", f"{pl(len(wraps), 'feature')} cross{'es' if len(wraps) == 1 else ''} the origin "
                       f"({', '.join(esc(f['name']) for f in wraps[:5])}); coordinates are written start–end through bp 1.")

    # ---- restriction
    mapping = cut_map(seq, circ)
    rows = enzyme_rows(rj, mapping)
    singles = [r for r in rows if r["cuts"] == 1]
    labels = single_cutter_labels(rj, mapping)
    mcs = [f for f in rj["features"] if re.search(r"\bMCS\b|multiple cloning|polylinker", f["name"] + " " + f["qualifiers"].get("note", ""), re.I)]
    for f in mcs[:1]:
        inside = [r["enzyme"] for r in singles if in_span(f, (r["positions"][0] - 1) % L, L)]
        if inside:
            R.flag("info", f"{pl(len(inside), 'single cutter')} inside <b>{esc(f['name'])}</b>: {', '.join(inside[:14])}{'…' if len(inside) > 14 else ''}.")
    R.tile(f"{commas(L)} bp", "length")
    R.tile(("circular" if circ else "linear"), "topology")
    R.tile(f"{gc:.1f}%", "GC content")
    R.tile(len(rj["features"]), "features")
    R.tile(len(singles), "single cutters (commercial)")
    R.r["summary"] = {"length": L, "circular": circ, "features": len(rj["features"]), "single_cutters": len(singles)}

    # ---- map
    R.section("map", "Construct", "Map", f"{esc(rj['name'])} · {commas(L)} bp · {'circular' if circ else 'linear'} · read from {esc(label)}.")
    orf0 = match_orfs(rj, find_orfs(seq, circ, 11, 75, "atg"))
    feat_types = {}
    for f in rj["features"]:
        feat_types[f["type"]] = feat_types.get(f["type"], 0) + 1
    R.fig("map", "Plasmid map" if circ else "Sequence map", "map",
          {"record_ref": True, "cutters": labels, "orfs": orf0[:60]},
          how=("Features are arrows in the direction they are read, coloured by type and labelled from /label, /gene, "
               "/product or /note. Enzymes that cut once (Biopython <code>Bio.Restriction</code>, commercially available "
               "enzymes only) are labelled around the outside at their cut site; isoschizomers that cut at the same "
               "position share one label. Hover for details; click a feature to select it and see its sequence."),
          yours=(f"{pl(len(rj['features']), 'feature')} ({', '.join(f'{v} {k}' for k, v in sorted(feat_types.items(), key=lambda x: -x[1])[:6])}) "
                 f"and {pl(len(singles), 'single cutter')} ({sum(1 for l in labels if l['six'] and l['neb'])} "
                 "positions with a 6+ bp site sold by NEB are shown by default)."),
          wide=True, sub="hover · click a feature · filter enzyme labels")

    # ---- features table
    R.section("features", "Annotation", "Features", "Coordinates are 1-based and inclusive, as in GenBank; a feature "
              "that crosses the origin of a circular plasmid is written start–end through bp 1.")
    frows = []
    for f in rj["features"]:
        chk = f.get("check")
        status = ""
        if chk:
            status = ("CHECK: " + "; ".join(chk["issues"])) if chk["issues"] else \
                ("OK" + (" · matches /translation" if chk["match"] else " · no /translation in file") +
                 (f" · {chk['aa']} aa" if chk.get("aa") else ""))
        frows.append([f["name"], f["type"], fmt_span(f), {1: "+", -1: "−"}.get(f["strand"], "·"), f["length"], status])
    R.table("features", "Features", ["Name", "Type", "Start–end (1-based)", "Strand", "Length (bp)", "CDS check"], frows,
            how=("For every CDS the DNA is translated (NCBI table from /transl_table, default 11) and compared with the "
                 "/translation qualifier; a mismatch, an internal stop, a missing start or stop, or a length that is not a "
                 "multiple of 3 is flagged CHECK. A valid alternative start codon (GTG, TTG) is read as Met, as ribosomes do."),
            yours=(f"{pl(len(cds_rows), 'CDS', 'CDS')} checked: {n_ok} clean, {n_bad} flagged." if cds_rows else
                   "No CDS features to check — a file without annotations (FASTA, plain sequence) has only the sequence."))

    # ---- restriction table
    R.section("restriction", "Bio.Restriction", "Restriction analysis",
              "Every commercially available enzyme in REBASE (Biopython <code>CommOnly</code>) searched with "
              f"<code>Analysis(..., linear={not circ})</code>. <b>Cut positions are the 1-based position of the first base "
              "after the cut on the top strand</b> (Biopython's convention). Left out: AbaSI, FspEI, LpnPI, MspJI and SgeI, which "
              "cut only methylated or hydroxymethylated DNA and so never cut a plasmid prep.")
    R.widget("enzymes", "enzymes", "Enzymes that cut", {"rows": rows, "noncutters": sorted([k for k, v in mapping.items() if not v], key=str.lower),
                                                         "n_searched": len(mapping)},
             how=("Filter to single or double cutters, to enzymes that cut once outside the selected feature "
                  "(<code>Analysis.do_not_cut</code> + one site), only inside it (<code>only_between</code>), the nearest "
                  "single cutters on each side of it, or enzymes that do not cut at all. Methylation sensitivity "
                  "(Dam, Dcm, CpG) is not modelled — check the supplier's page for enzymes near GATC or CCWGG."),
             yours=(f"{sum(1 for v in mapping.values() if v)} of {len(mapping)} commercial enzymes cut this sequence; "
                    f"{len(singles)} cut exactly once and {sum(1 for r in rows if r['cuts'] == 2)} twice."))

    # ---- digest
    R.section("digest", "Virtual digest", "Digest and gel", "Pick up to three enzymes per lane; fragments respect "
              f"topology ({'n cuts give n fragments on a circle' if circ else 'n cuts give n + 1 fragments on a line'}).")
    sug = []
    try:
        sug = suggest_digests(rj, None, 1.0, top=3)
    except Exception:  # noqa: BLE001
        sug = []
    lanes = [[]]
    single_names = {r["enzyme"] for r in singles}
    first_single = next((n for n in COMMON if n in single_names), None) or \
        next((r["enzyme"] for r in singles if enzyme_info(enzyme(r["enzyme"]))["neb"] and r["size"] >= 6), None)
    if first_single:
        lanes.append([first_single])
    if sug:
        lanes.append(sug[0]["enzymes"])
    dg = digest_api({"record": rj, "lanes": lanes, "agarose": 1.0})
    R.widget("digest", "digest", "Virtual agarose gel", {"initial": dg, "suggestions": sug},
             how=("Band positions follow a log-linear size–migration relation inside the range each agarose percentage "
                  "resolves (0.7 %: ~0.8–12 kb, 1 %: ~0.5–10 kb, 2 %: ~0.1–2 kb), compressed above it. Band brightness scales "
                  "with the mass of DNA in it (equimolar fragments: bigger is brighter). Fragments within ~5 % of each other's "
                  "size co-migrate as one band; fragments under ~100 bp are faint or run off. Uncut circular plasmid is "
                  "drawn without a size: supercoiled DNA runs faster than its length, and nicked/linear forms run elsewhere."),
             yours=dg["yours"], sub="add lanes · change ladder and agarose · compare with the construct minus a feature")

    # ---- ORFs
    R.section("orfs", "Six-frame scan", "Open reading frames", "Start codon to the next in-frame stop in all six frames"
              + (", scanning across the origin" if circ else "") + ".")
    R.widget("orfs", "orfs", "ORFs", {"initial": {"orfs": orf0[:300], "n": len(orf0), "table": 11, "min_aa": 75, "starts": "atg",
                                                  "yours": orf_yours(rj, orf0, 11, 75, "atg")},
                                      "tables": codon_tables_api()["tables"]},
             how=("Each frame is read codon by codon with the chosen NCBI genetic code (<code>Bio.Data.CodonTable</code>); an ORF "
                  "runs from the first start codon after a stop to the next stop. Circular sequences are scanned as "
                  "sequence + sequence so ORFs crossing the origin are found once. Proteins come from "
                  "<code>Seq.translate(table=…)</code>. “= name” means the ORF has the same start and stop as an annotated CDS."),
             yours=orf_yours(rj, orf0, 11, 75, "atg"))

    R.method("Reading", f"{esc(label)} read with Biopython {__import__('Bio').__version__} <code>SeqIO</code>"
                        f"{' (binary mode)' if label in ('SnapGene .dna', 'DNA Strider / Serial Cloner .xdna', 'Gene Construction Kit .gck') else ''}; "
                        "topology from <code>annotations['topology']</code>.")
    R.method("CDS checks", "Each CDS was extracted with its location (compound locations across the origin included), "
                           "translated with <code>Seq.translate</code> using its /transl_table (default 11) and compared residue by "
                           "residue with /translation.")
    R.method("Restriction analysis", "<code>Bio.Restriction.Analysis</code> with the REBASE commercially available set "
                                     f"(<code>CommOnly</code> without the five modification-dependent enzymes, {len(mapping)} enzymes), "
                                     f"linear={not circ}. Cut positions: 1-based "
                                     "first base after the top-strand cut. REBASE: Roberts et al., Nucleic Acids Res 2023.")
    R.method("Virtual gel", "Fragment sizes from the cut positions; migration drawn log-linear in size within the resolvable "
                            "range for the agarose percentage — an illustration, not a calibrated model.")
    R.method("ORFs", "Six-frame scan with <code>Bio.Data.CodonTable</code> (NCBI genetic codes) and <code>Seq.translate</code>.")
    return R.done()


def translate_feature_api(args):
    rj = args.get("record")
    if not rj:
        raise UserFacingError("Load a construct first.")
    f = find_feature(rj, args.get("feature"))
    if f is None:
        raise UserFacingError(f"No feature called “{args.get('feature')}”. Features: "
                              + ", ".join(x["name"] for x in rj["features"][:30]))
    nt = extract_feature(rj["seq"], f)
    chk = check_cds(f, rj["seq"]) if f["type"] == "CDS" else None
    if chk is None:
        table = int(args.get("table") or 11)
        usable = nt[:len(nt) // 3 * 3]
        prot = str(Seq(usable).translate(table=table))
    else:
        prot = chk["translation"]
    return {"feature": f["name"], "type": f["type"], "span": fmt_span(f), "strand": f["strand"], "length_bp": len(nt),
            "dna": nt, "protein": prot, "check": chk}


# ================================================================ exports
def export_api(args):
    rj = args.get("record")
    fmt = (args.get("format") or "genbank").lower()
    rec = record_from_json(rj)
    rec.name = re.sub(r"\W+", "_", rec.name)[:16] or "construct"
    if fmt == "fasta":
        return {"text": f">{rec.id} {rj.get('name', '')} {len(rec)} bp {'circular' if rj.get('circular') else 'linear'}\n"
                        + "\n".join(str(rec.seq)[i:i + 70] for i in range(0, len(rec), 70)) + "\n",
                "filename": re.sub(r"\W+", "_", rj.get("name") or "construct") + ".fasta"}
    for f in rec.features:
        f.qualifiers.pop("translation", None) if f.type != "CDS" else None
    buf = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SeqIO.write(rec, buf, "genbank")
    return {"text": buf.getvalue(), "filename": re.sub(r"\W+", "_", rj.get("name") or "construct") + ".gb"}


def read_reference_api(args):
    """A reference sequence for Sanger / PCR checks, from a file or pasted text."""
    if args.get("file"):
        rec, label, notes = parse_sequence_file(args["file"], None if args.get("circular") in (None, "", "auto") else bool(args.get("circular")))
    else:
        rec, notes = record_from_text(args.get("text") or "", bool(args.get("circular")), args.get("name") or "Pasted reference")
    rj = record_json(rec)
    for f in rj["features"]:
        if f["type"] == "CDS":
            f["translation"] = check_cds(f, rj["seq"])["translation"]
    return {"record": rj, "notes": notes}


# ================================================================ sequence tools
def booth_min_rotation(s: str) -> int:
    """Index of the lexicographically smallest rotation (Booth 1980), O(n)."""
    S = s + s
    f = [-1] * len(S)
    k = 0
    for j in range(1, len(S)):
        sj = S[j]
        i = f[j - k - 1]
        while i != -1 and sj != S[k + i + 1]:
            if sj < S[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != S[k + i + 1]:
            if sj < S[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return k % len(s)


def canonical_circular(s: str) -> tuple[str, int, int]:
    """(canonical sequence, strand, rotation): the smallest rotation over both strands."""
    r1 = booth_min_rotation(s)
    c1 = s[r1:] + s[:r1]
    rc = revcomp(s)
    r2 = booth_min_rotation(rc)
    c2 = rc[r2:] + rc[:r2]
    return (c1, 1, r1) if c1 <= c2 else (c2, -1, r2)


def seq_input(args, key="text", fkey="file", circ_key="circular"):
    if args.get(fkey):
        rec, _, notes = parse_sequence_file(args[fkey], None if args.get(circ_key) in (None, "", "auto") else bool(args.get(circ_key)))
        return str(rec.seq), rec.annotations.get("topology") == "circular", args[fkey].get("name") or rec.id, notes
    if args.get(key):
        rec, notes = record_from_text(args[key], bool(args.get(circ_key)), "Sequence")
        return str(rec.seq), bool(args.get(circ_key)), rec.name, notes
    if args.get("record"):
        rj = args["record"]
        return rj["seq"], rj["circular"], rj["name"], []
    raise UserFacingError("Paste a sequence or drop a file.")


def motif_search(seq: str, motif: str, circular: bool) -> list:
    m, _ = normalize_dna(motif, "The motif")
    L, n = len(seq), len(m)
    if n > L:
        return []
    ext = seq + seq[:n - 1] if circular else seq
    hits = []
    fw = nt_search(ext, m)[1:]
    for p in fw:
        if p < L:
            hits.append({"strand": 1, "start0": p, "end0": (p + n - 1) % L + 1, "match": ext[p:p + n]})
    if revcomp(m) != m:
        rc = revcomp(ext)
        E = len(ext)
        for p in nt_search(rc, m)[1:]:
            f0 = E - p - n            # forward-strand start of the reverse-strand hit
            if 0 <= f0 < L:
                hits.append({"strand": -1, "start0": f0, "end0": (f0 + n - 1) % L + 1, "match": revcomp(ext[f0:f0 + n])})
    hits.sort(key=lambda h: (h["start0"], -h["strand"]))
    return hits


def seqtools(args: dict) -> dict:
    seq, circ, name, notes = seq_input(args)
    L = len(seq)
    table = int(args.get("table") or 1)
    R = Result("seqtools", name)
    R.r["record"] = {"name": name, "length": L, "circular": circ, "seq": seq, "features": []}
    for lvl, t in notes:
        R.flag(lvl, t)
    R.tile(f"{commas(L)} bp", "length")
    R.tile("circular" if circ else "linear", "topology")
    R.tile(f"{100 * gc_fraction(seq, ambiguous='ignore'):.1f}%", "GC")
    sg = seguid(seq)
    can, cstrand, crot = canonical_circular(seq) if circ else (seq, 1, 0)
    csg = seguid(can) if circ else None
    R.tile(sg, "SEGUID (as given)")
    R.section("rc", "Strands", "Reverse complement", "")
    R.fig("rc", "Reverse complement", "mono", {"text": revcomp(seq), "wrap": 60},
          how="Complement of each base (IUPAC codes included), read 5′→3′.", yours=f"{commas(L)} bp.")
    R.section("tr", "Translation", f"Translation · NCBI table {table}", "")
    frames = []
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        for fr in range(3):
            sub = s[fr:]
            sub = sub[:len(sub) // 3 * 3]
            frames.append({"frame": f"{'+' if strand == 1 else '−'}{fr + 1}", "protein": str(Seq(sub).translate(table=table))})
    R.fig("translate", "Six-frame translation", "frames", {"frames": frames},
          how=f"<code>Seq.translate(table={table})</code> of each frame; * is a stop codon.",
          yours=f"Frame +1 reads {frames[0]['protein'][:40]}{'…' if len(frames[0]['protein']) > 40 else ''} "
                f"({frames[0]['protein'].count('*')} stops).")
    motif = (args.get("motif") or "").strip()
    if motif:
        hits = motif_search(seq, motif, circ)
        R.section("motif", "IUPAC search", f"Motif {esc(motif.upper())}", "")
        R.table("motif", "Motif hits (both strands)", ["Strand", "Start (1-based)", "End", "Matched sequence"],
                [["+" if h["strand"] == 1 else "−", h["start0"] + 1, h["end0"], h["match"]] for h in hits],
                how=("<code>nt_search</code> on the sequence and on its reverse complement; reverse-strand hits are converted "
                     "to forward-strand coordinates. IUPAC codes in the motif (N, R, Y, W…) match any of their bases."
                     + (" The circular sequence is searched across the origin." if circ else "")),
                yours=f"{pl(len(hits), 'hit')} for {esc(motif.upper())}" + (" (the motif is its own reverse complement, so each site is listed once)." if revcomp(motif.upper()) == motif.upper() else "."))
    R.section("seguid", "Checksum", "Same sequence as another file?", "")
    cmp_rows = [["This sequence", name, commas(L), sg, csg or "—"]]
    yours = f"SEGUID {sg}."
    other = None
    if args.get("other_text") or args.get("other_file") or args.get("other_record"):
        oa = {"text": args.get("other_text"), "file": args.get("other_file"), "record": args.get("other_record"),
              "circular": args.get("other_circular", circ)}
        oseq, ocirc, oname, _ = seq_input(oa)
        osg = seguid(oseq)
        ocan = canonical_circular(oseq) if (circ or ocirc) else (oseq, 1, 0)
        ocsg = seguid(ocan[0]) if (circ or ocirc) else None
        cmp_rows.append(["Other sequence", oname, commas(len(oseq)), osg, ocsg or "—"])
        if osg == sg:
            yours = "<b>Identical</b>: the two sequences have the same SEGUID, base for base, from the same starting point."
            R.flag("ok", "The two sequences are identical.")
        elif (circ or ocirc) and len(oseq) == L and csg and seguid(canonical_circular(seq)[0]) == ocsg:
            # the same circle: find the rotation and strand that maps one onto the other
            ss = seq + seq
            off = ss.find(oseq)
            strand = 1
            if off < 0:
                off = (seq + seq).find(revcomp(oseq))
                strand = -1
            yours = ("<b>The same circular molecule</b>, written from a different starting point"
                     + (" and on the other strand" if strand == -1 else "")
                     + ((f": the other sequence{'’s reverse complement' if strand == -1 else ''} starts at bp "
                         f"{commas(off + 1)} of this one") if off >= 0 else "") + ".")
            R.flag("ok", "Same circular plasmid, rotated" + (" and reverse-complemented" if strand == -1 else "") + ".")
        else:
            diff = len(oseq) - L
            yours = (f"<b>Different sequences</b> ({commas(L)} vs {commas(len(oseq))} bp"
                     + (f", {abs(diff)} bp {'longer' if diff > 0 else 'shorter'}" if diff else ", same length") + ").")
            R.flag("info", "The two sequences differ. Align them (Sanger tool, or any aligner) to see where.")
        other = True
    R.table("seguid", "Checksums", ["", "Name", "Length", "SEGUID (as given)", "Circular SEGUID (rotation- and strand-free)"],
            cmp_rows,
            how=("SEGUID is a SHA-1 of the upper-case sequence (Biopython <code>CheckSum.seguid</code>, Babnigg & Giometti 2006): "
                 "identical sequences share it. A plasmid can be written from any starting point and on either strand, so for "
                 "circular sequences the app also hashes a canonical form — the lexicographically smallest rotation over both "
                 "strands (Booth's algorithm). This is this app's own convention, not the published cSEGUID (SEGUID v2) standard."),
            yours=yours if other else yours + " Add a second sequence to compare.")
    R.method("Tools", "Reverse complement and translation with Biopython <code>Seq</code>; motif search with "
                      "<code>Bio.SeqUtils.nt_search</code> on both strands; checksums with <code>Bio.SeqUtils.CheckSum.seguid</code>.")
    return R.done()
