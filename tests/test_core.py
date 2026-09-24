"""Core tests: numbers checked against known values or against Biopython called directly.

    .venv/bin/python -m pytest tests -q
"""
import base64
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from Bio import SeqIO  # noqa: E402
from Bio.Restriction import EcoRI  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqUtils import MeltingTemp as mt  # noqa: E402

from server import core_primers, core_protein, core_sanger, core_seq  # noqa: E402
from server.core import UserFacingError, api  # noqa: E402

EX = ROOT / "examples"


def filearg(name, data=None):
    data = (EX / name).read_bytes() if data is None else data
    return {"name": name, "b64": base64.b64encode(data).decode()}


fileArg = filearg


def rec_json(seq, circular, features=()):
    return {"name": "t", "id": "t", "length": len(seq), "circular": circular, "seq": seq, "features": list(features)}


# ---------------------------------------------------------------- primers / Tm
def test_tm_nn_passes_the_stated_conditions():
    primer = "AGCGGATAACAATTTCACACAGGA"            # M13 reverse primer
    b = core_primers.buffer_from({"preset": "taq"})
    expect = mt.Tm_NN(primer, nn_table=mt.DNA_NN3, saltcorr=7, Na=0, K=50, Tris=10, Mg=1.5, dNTPs=0.8, dnac1=250, dnac2=0)
    assert core_primers.tm_nn(primer, b) == pytest.approx(expect, abs=1e-9)
    # and through the whole tool, rounded as displayed
    r = api("primers", {"text": f"M13rev {primer}", "preset": "taq"})
    assert r["primers"][0]["tm_nn"] == pytest.approx(expect, abs=1e-9)
    # DMSO goes through chem_correction (0.75 °C per %)
    b5 = core_primers.buffer_from({"preset": "taq", "buffer": {"DMSO": 5}})
    assert core_primers.tm_nn(primer, b5) == pytest.approx(mt.chem_correction(expect, DMSO=5), abs=1e-9)
    assert b5["label"].startswith("Custom")
    # the high-fidelity preset: 2 mM Mg, 500 nM primer, no monovalent salt
    bh = core_primers.buffer_from({"preset": "hifi"})
    eh = mt.Tm_NN(primer, nn_table=mt.DNA_NN3, saltcorr=7, Na=0, K=0, Tris=0, Mg=2.0, dNTPs=0.8, dnac1=500, dnac2=0)
    assert core_primers.tm_nn(primer, bh) == pytest.approx(eh, abs=1e-9)
    assert "proprietary" in bh["why"] and "vendor calculator" in bh["why"]


def test_primer_checks_and_heuristics():
    r = core_primers.primer_row({"name": "p", "raw": "GGGGCCCCGGGGCCCCGCGC"}, core_primers.buffer_from({}))
    assert r["clamp_level"] == "warn"                    # 5 G/C in the last five bases
    assert any("run of 4" in x for x in r["repeats"])
    sd = core_primers.comp_scan("GAATTCGAATTC", "GAATTCGAATTC")   # palindromic: fully self-complementary
    assert sd["longest"]["run"] == 12 and sd["three"]["run"] >= 3
    hp = core_primers.hairpin_scan("AAAAGGGCCCTTTTTTTGGGCCCAAAA")
    assert hp["stem"] >= 6


def test_pcr_across_the_origin_of_a_circular_template():
    import random
    random.seed(1)
    seq = "".join(random.choice("ACGT") for _ in range(3000))
    fwd = seq[2900:2922]                                  # forward primer near the end
    rev = core_seq.revcomp(seq[80:102])                   # reverse primer near the start
    res = core_primers.pcr([{"name": "F", "seq": fwd}, {"name": "R", "seq": rev}], rec_json(seq, True))
    main = [p for p in res["products"] if p["perfect"]][0]
    assert main["wraps"] and main["size"] == (3000 - 2900) + 102


# ---------------------------------------------------------------- restriction & digest
def test_ecori_linear_and_circular():
    seq = "A" * 100 + "GAATTC" + "T" * 200 + "GAATTC" + "C" * 94       # 406 bp
    lin = rec_json(seq, False)
    m = core_seq.cut_map(seq, False, ["EcoRI"])
    assert m["EcoRI"] == [102, 308]                      # 1-based first base after G^AATTC
    assert m["EcoRI"] == EcoRI.search(Seq(seq), linear=True)
    d = core_seq.digest(lin, ["EcoRI"])
    assert sorted(f["size"] for f in d["fragments"]) == sorted([101, 206, 99])     # n cuts → n + 1 fragments
    circ = rec_json(seq, True)
    dc = core_seq.digest(circ, ["EcoRI"])
    assert sorted(f["size"] for f in dc["fragments"]) == [200, 206]                 # n cuts → n fragments
    # a site across the origin is found only when the sequence is circular
    wrap = "AATTC" + "A" * 50 + "G"
    assert core_seq.cut_map(wrap, True, ["EcoRI"])["EcoRI"] == [1]
    assert core_seq.cut_map(wrap, False, ["EcoRI"])["EcoRI"] == []
    # uncut circular plasmid: no size drawn
    u = core_seq.digest(rec_json("A" * 500, True), ["EcoRI"])
    assert u["uncut_circular"] and u["bands"][0]["size"] is None


def test_comigration_and_filters_on_pbad30():
    r = api("construct", {"file": fileArg("pBAD30.gb")})
    rj = r["record"]
    bands = core_seq.bands_for([{"size": 1000}, {"size": 1030}, {"size": 3000}], 1.0)
    assert [b["n"] for b in bands] == [1, 2]            # 1,000 and 1,030 bp co-migrate
    araC = core_seq.find_feature(rj, "araC")
    out = core_seq.filter_enzymes(rj, core_seq.cut_map(rj["seq"], True), "outside", "araC")
    for row in out["rows"]:
        assert row["cuts"] == 1
        p0 = row["positions"][0] - 1
        assert not (araC["start0"] < p0 < araC["end0"])
    assert "do_not_cut" in out["note"]


# ---------------------------------------------------------------- reading & honest CDS checks
@pytest.mark.parametrize("name,length,topo", [("pBAD30.gb", 4923, True), ("addgene-plasmid-39296-sequence-49545.gbk", 3938, True),
                                              ("pFA-KanMX4.dna", 3941, True)])
def test_examples_read(name, length, topo):
    r = api("construct", {"file": fileArg(name)})
    assert r["record"]["length"] == length and r["record"]["circular"] is topo


def test_cds_check_finds_the_pbad30_annotation_slip_and_passes_good_ones():
    r = api("construct", {"file": fileArg("pBAD30.gb")})
    apr = core_seq.find_feature(r["record"], "AP(R)")
    assert apr["check"]["issues"] and "runs past its stop codon" in apr["check"]["issues"][0]
    g = api("construct", {"file": fileArg("addgene-plasmid-39296-sequence-49545.gbk")})
    for n in ("KanR", "AmpR"):
        f = core_seq.find_feature(g["record"], n)
        assert f["check"]["issues"] == [] and f["check"]["match"] is True
    # a wrong /translation is caught
    f = dict(core_seq.find_feature(g["record"], "KanR"))
    f["qualifiers"] = dict(f["qualifiers"], translation="MSIQ" + f["qualifiers"]["translation"][4:])
    chk = core_seq.check_cds(f, g["record"]["seq"])
    assert chk["match"] is False and any("does not match" in x for x in chk["issues"])


def test_origin_wrapping_feature():
    g = api("construct", {"file": fileArg("addgene-plasmid-39296-sequence-49545.gbk")})
    sp6 = core_seq.find_feature(g["record"], "SP6 promoter")
    assert sp6["wraps"] and sp6["start0"] == 3921 and sp6["end0"] == 2


def test_bad_inputs_have_clear_messages():
    with pytest.raises(UserFacingError, match="Sanger trace"):
        api("construct", {"file": fileArg("3730.ab1")})
    with pytest.raises(UserFacingError, match="empty"):
        api("construct", {"file": {"name": "empty.gb", "b64": ""}})
    with pytest.raises(UserFacingError, match="protein"):
        api("construct", {"text": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", "circular": False})
    with pytest.raises(UserFacingError, match="not a Sanger trace"):
        api("sanger", {"files": [fileArg("pBAD30.gb")], "reference_text": "ACGT" * 50})
    rna = api("construct", {"text": "AUGGCUAGCAAAGGAGAAGAACUUUUCACUGGAGUUGUCCCAAUUCUUGUUGAAUUAGAUGG", "circular": False})
    assert "RNA" in " ".join(f["text"] for f in rna["flags"]) and "U" not in rna["record"]["seq"]
    soft = api("construct", {"text": "acgtacgtacGAATTCacgtacgtacgtNNNNNNNNNNacgtacgtacgtacgtacgt", "circular": False})
    txt = " ".join(f["text"] for f in soft["flags"])
    assert "lower-case" in txt and "<b>N</b>" in txt
    assert soft["record"]["seq"].isupper()


# ---------------------------------------------------------------- ORFs
def test_orf_across_the_origin():
    body = "ATG" + "GCT" * 99 + "TAA"                    # 100 codons + stop = 303 bp
    spacer = "C" * 200
    seq = body[150:] + spacer + body[:150]               # the ORF starts near the end and wraps
    orfs = core_seq.find_orfs(seq, True, 11, 75, "atg")
    hit = [o for o in orfs if o["strand"] == 1 and o["aa"] == 100]
    assert hit and hit[0]["wraps"] and hit[0]["start0"] == len(seq) - 150 and hit[0]["end0"] == 153
    assert not [o for o in core_seq.find_orfs(seq, False, 11, 75, "atg") if o["aa"] == 100]


# ---------------------------------------------------------------- Sanger
def test_mott_trim_matches_biopython_abi_trim():
    tr = core_sanger.read_trace(fileArg("3730.ab1"))
    s, e = core_sanger.mott_trim(tr["quals"])
    bio = str(SeqIO.read(EX / "3730.ab1", "abi-trim").seq)
    assert tr["calls"][s:e] == bio.upper()


def test_planted_substitution_and_deletion_are_found_with_consequences():
    ab1 = fileArg("3730.ab1")
    d = core_sanger.build_demo_reference(ab1)
    rj = core_sanger.demo_record(d)
    r = api("sanger", {"files": [ab1], "reference": rj})
    assert r["summary"]["confident"] == 2
    read = core_sanger.analyse_read(core_sanger.read_trace(ab1), rj)
    subs = [x for x in read["diffs"] if x["type"] == "substitution"]
    dels = [x for x in read["diffs"] if x["type"] == "deletion"]
    assert len(subs) == 1 and len(dels) == 1 and len(read["diffs"]) == 2
    assert subs[0]["ref_pos1"] == d["substitution"]["ref_pos1"]
    assert subs[0]["ref"] == d["substitution"]["ref"] and subs[0]["read"] == d["substitution"]["read"]
    assert subs[0]["confidence"] == "confident"
    eff = subs[0]["effects"][0]
    codon = (d["substitution"]["ref_pos1"] - 1 - d["atg"]) // 3 + 1
    assert eff["kind"] == "missense" and eff["codon"] == codon
    assert str(Seq(rj["seq"][d["atg"] + 3 * (codon - 1): d["atg"] + 3 * codon]).translate()) == eff["short"][0]
    assert dels[0]["ref_pos1"] == d["extra_base"]["ref_pos1"] and dels[0]["ref"] == d["extra_base"]["base"]
    assert dels[0]["effects"][0]["kind"] == "frameshift"
    # the saved demo files are the same reference
    fa = "".join(l.strip() for l in (EX / "demo_reference.fasta").read_text().splitlines()[1:])
    assert fa == d["seq"]


def test_consequence_on_a_minus_strand_cds():
    # a tiny minus-strand CDS: ATG GCT TGG TAA on the bottom strand
    cds = "ATGGCTTGGTAA"
    seq = "CCCC" + core_seq.revcomp(cds) + "CCCC"
    fj = {"name": "x", "type": "CDS", "strand": -1, "parts": [[4, 16]], "qualifiers": {}}
    # top-strand position of the W codon's middle G (cds index 7) is 4 + (11 - 7) = 8; top base C → read T means cds G → A: TGG→TAG
    eff = core_sanger.codon_change(fj, seq, 8, "T")
    assert eff["kind"] == "nonsense" and eff["codon"] == 3


def test_310_trace_without_quality_is_reported_not_trusted():
    tr = core_sanger.read_trace(fileArg("310.ab1"))
    assert tr["has_q"] is False
    ref = "".join(c if c in "ACGT" else "A" for c in tr["calls"][100:500])
    rj = rec_json(ref, False)
    read = core_sanger.analyse_read(tr, rj)
    assert read["aligned"]
    assert all(x["confidence"] == "unscored" for x in read["diffs"])


# ---------------------------------------------------------------- protein
LYSOZYME = ("KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNC"
            "AKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")


def test_lysozyme_extinction_coefficient():
    r = api("protein", {"text": LYSOZYME, "mode": "protein"})
    eps = r["protein"]["ext_reduced"]
    assert 36000 <= eps <= 38000 and eps == 6 * 5500 + 3 * 1490
    assert r["protein"]["ext_cystines"] == eps + 4 * 125


def test_aa_percent_shim_both_biopython_behaviours():
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(LYSOZYME)
    new = core_protein.compat_aa_percent(pa)             # 1.87+: cached property, percent
    assert sum(new.values()) == pytest.approx(100, abs=0.01)

    class Old185:                                        # Biopython 1.85: method returning fractions,
        amino_acids_percent = None                       # plain attribute that is None until it is called

        def __init__(self, seq):
            self.sequence = seq

        def get_amino_acids_percent(self):
            return {a: self.sequence.count(a) / len(self.sequence) for a in core_protein.STD}

    old = core_protein.compat_aa_percent(Old185(LYSOZYME))
    for a in core_protein.STD:
        assert old[a] == pytest.approx(new[a], abs=1e-9)


def test_tags_and_scars():
    r = api("protein", {"text": "MGSSHHHHHHSSGLVPRGSHMENLYFQG" + LYSOZYME, "mode": "protein"})
    rows = [it for s in r["sections"] for it in s["items"] if it["id"] == "tags"][0]["rows"]
    names = [x[0] for x in rows]
    assert "His-tag" in names and "TEV site" in names and "Thrombin site" in names
    tev = [x for x in rows if x[0] == "TEV site"][0]
    assert "starting GKVFGR" in tev[4]


# ---------------------------------------------------------------- sequence tools
def test_seguid_recognises_a_rotated_reverse_complemented_plasmid():
    seq = str(SeqIO.read(EX / "pBAD30.gb", "genbank").seq)
    rot = core_seq.revcomp(seq[1234:] + seq[:1234])
    r = api("seqtools", {"text": seq, "circular": True, "other_text": rot, "other_circular": True})
    assert any("Same circular plasmid" in f["text"] for f in r["flags"])
    r2 = api("seqtools", {"text": seq, "circular": True, "other_text": seq[:-1] + "A", "other_circular": True})
    assert not any("Same circular" in f["text"] for f in r2["flags"])


def test_motif_positions_on_both_strands():
    seq = "AAAGGTCTCAAAAAAGAGACCAAA"                     # BsaI site forward at 4, reverse at 16
    hits = core_seq.motif_search(seq, "GGTCTC", False)
    assert [(h["strand"], h["start0"] + 1) for h in hits] == [(1, 4), (-1, 16)]


def test_genbank_export_round_trip(tmp_path):
    g = api("construct", {"file": fileArg("addgene-plasmid-39296-sequence-49545.gbk")})
    txt = api("export", {"record": g["record"], "format": "genbank"})["text"]
    p = tmp_path / "x.gb"
    p.write_text(txt)
    rec = SeqIO.read(p, "genbank")
    assert str(rec.seq) == g["record"]["seq"] and rec.annotations["topology"] == "circular"
    assert rec.annotations["molecule_type"] == "DNA"


def test_primer_names_that_look_like_bases():
    """'R ACGT…' and 'KanR TTAG…' are names: every letter of R and KanR is an IUPAC code, and read as sequence they
    glued degenerate 'bases' onto the 5′ end (found by randomised testing: the product came out 1–4 bp long)."""
    from server.core_primers import parse_primers
    got = [(p["name"], p["raw"]) for p in parse_primers(
        "F ACGTACGTACGTAGCTAG\nR TTGCATGCATGCAAGCTT\nKanR TTAGAAAAACTCATCGAGCATC\nACGTACGTAC GTACGTACGT")]
    assert got == [("F", "ACGTACGTACGTAGCTAG"), ("R", "TTGCATGCATGCAAGCTT"),
                   ("KanR", "TTAGAAAAACTCATCGAGCATC"), ("P4", "ACGTACGTACGTACGTACGT")]


def test_modification_dependent_enzymes_left_out():
    """AbaSI, MspJI… cut only methylated DNA; as plain patterns they 'cut' at almost every C of a plasmid."""
    from server.core_seq import _batch
    names = {str(e) for e in _batch()}
    assert not names & {"AbaSI", "FspEI", "LpnPI", "MspJI", "SgeI"} and "EcoRI" in names


def test_non_ascii_sequence_is_a_plain_error():
    """A symbol pasted from Word/PDF into a FASTA crashed with UnicodeDecodeError (found by randomised testing)."""
    import pytest
    from server import core
    from server.common import UserFacingError
    with pytest.raises(UserFacingError, match="not DNA letters"):
        core.api("construct", {"file": {"name": "x.fasta", "text": ">x\nACGTΩ≈çACGTACGTACGTACGTACGTACGTACGTACGT\n"}})


def test_pdf_report_for_every_tool(tmp_path):
    """Every tool's result must make a PDF. The Sanger one never could: a bare '<' in its methods text ('Q < 20') was
    read as a tag and broke ReportLab's markup (found by repeated HTTP testing). Also checks the date_toolname_report name."""
    import base64
    import os
    import re
    os.environ["CB_HOME"] = str(tmp_path)
    from server import core
    from server.common import report_filename, save_item
    from server.report import build_pdf
    ex = os.path.join(os.path.dirname(__file__), "..", "examples")
    f = lambda n: {"name": n, "b64": base64.b64encode(open(os.path.join(ex, n), "rb").read()).decode()}  # noqa: E731
    runs = {"construct": {"file": f("pBAD30.gb")}, "primers": {"text": "F ATGGGTAAGGAAAAGACTCACG\nR TTAGAAAAACTCATCGAGCATC"},
            "sanger": {"files": [f("3730.ab1")], "reference_file": f("demo_reference.gb"), "demo": True},
            "protein": {"text": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"}, "seqtools": {"text": "ACGTGAATTCAAGAATTCTTACGATCG", "motif": "GAATTC"}}
    for kind, args in runs.items():
        res = core.api(kind, args)
        iid = save_item(kind, args, res)
        out = tmp_path / f"{kind}.pdf"
        build_pdf(res, out)
        assert out.read_bytes()[:5] == b"%PDF-", kind
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}_CloneBench-[A-Za-z]+_report\.pdf", report_filename(iid, res)), kind


# ---------------------------------------------------------------- primer design (0.2.0)
def _kan():
    import base64
    import os
    from server import core
    ex = os.path.join(os.path.dirname(__file__), "..", "examples", "addgene-plasmid-39296-sequence-49545.gbk")
    return core.api("construct", {"file": {"name": "kan.gbk", "b64": base64.b64encode(open(ex, "rb").read()).decode()}})["record"]


def test_design_pcr_pair_is_on_the_template_and_amplifies_the_target():
    from server import core
    from server.core_primers import pcr
    rj = _kan()
    r = core.api("design", {"template": rj, "mode": "pcr", "feature": "KanR"})
    p = r["design"]["pairs"][0]
    s = rj["seq"].upper()
    assert s[p["fwd_start0"]:p["fwd_start0"] + len(p["fwd"])] == p["fwd"]
    prods = pcr([{"name": "F", "seq": p["fwd"]}, {"name": "R", "seq": p["rev"]}], rj)["products"]
    assert any(x["size"] == p["size"] and x["perfect"] for x in prods)
    t = r["design"]["target"]
    assert p["fwd_start0"] <= t["start0"] and p["rev_end0"] >= t["start0"] + t["length"]


def test_design_ncoi_merges_the_start_codon_and_reports_a_codon_change():
    """NcoI (CCATGG) in front of ATG AAA… must become …CC ATG GAA… (K2E, reported), never CCATGG + ATGAAA (an
    out-of-frame ATG ahead of the gene in a pET vector)."""
    from server import core
    from server.core_seq import record_from_text, record_json
    cds = "ATGAAAGCTTTCGGTACCGGAGCTAGCCTGGATGTTAAACCGGCATTTGGCTATCGTCGTGGCAAAGATCTGGCCGTAA"
    rj = record_json(record_from_text("GGGCCCTTTAAAGGGCCC" + cds + "CCCGGGTTTAAACCCGGG", False, "t")[0])
    rj["features"].append({"i": 0, "name": "myCDS", "type": "CDS", "strand": 1, "parts": [[18, 18 + len(cds)]], "start0": 18,
                           "end0": 18 + len(cds), "wraps": False, "length": len(cds), "qualifiers": {}})
    r = core.api("design", {"template": rj, "mode": "clone", "feature": "myCDS", "enzyme5": "NcoI", "enzyme3": "KpnI"})
    assert r["design"]["fwd"].startswith("gcgcccATGGAAGCT")
    assert any("K2E" in f["text"] for f in r["flags"])
    assert any("KpnI cuts inside the insert" in f["text"] for f in r["flags"])      # GGTACC is in the CDS
    r = core.api("design", {"template": rj, "mode": "clone", "feature": "myCDS", "enzyme5": "NdeI"})
    assert r["design"]["fwd"].startswith("gcgccatATGAAAGCT")


def test_design_sequencing_walk_covers_the_target():
    from server import core
    rj = _kan()
    r = core.api("design", {"template": rj, "mode": "seq", "feature": "kanMX", "read_len": 600})
    t, L = r["design"]["target"], rj["length"]
    covered = set()
    for rd in r["design"]["reads"]:
        off = ((rd["from0"] - t["start0"] + L // 2) % L) - L // 2
        covered |= set(range(off, off + rd["len"]))
    assert set(range(t["length"])) <= covered
