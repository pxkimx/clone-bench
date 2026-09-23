"""Quick end-to-end check of this installation, offline:  python -m server.selftest

Runs every tool on the bundled examples, saves each as a workspace item in a temporary CB_HOME, writes a PDF
report, and exercises the assistant's tools (without calling the API). Prints pass/fail; exit code 1 on failure.
"""
from __future__ import annotations

import base64
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="cb_selftest_"))
    os.environ["CB_HOME"] = str(tmp / "home")
    from . import agent, core
    from .common import ROOT
    from .common import load_item, save_item
    from .report import build_pdf

    ex = ROOT / "examples"
    f = lambda n: {"name": n, "b64": base64.b64encode((ex / n).read_bytes()).decode()}  # noqa: E731
    print(f"Clone Bench self-test · Biopython {core.versions()['biopython']} · python {core.versions()['python']}")
    ok = True
    items = {}

    def check(label, name, args, expect=None):
        nonlocal ok
        t = time.time()
        try:
            r = core.api(name, args)
            if expect:
                expect(r)
            if name in ("construct", "primers", "sanger", "protein", "seqtools"):
                items[label] = save_item(name, args, r)
            tiles = " · ".join(f"{x['value']} {x['label']}" for x in r.get("tiles", [])[:4])
            print(f"  ok   {label:<34} {time.time() - t:5.2f}s  {tiles}")
            return r
        except Exception as e:  # noqa: BLE001
            ok = False
            print(f"  FAIL {label:<34} {type(e).__name__}: {e}")
            traceback.print_exc()
            return None

    def expect_error(label, name, args, words):
        nonlocal ok
        try:
            core.api(name, args)
            ok = False
            print(f"  FAIL {label:<34} no error raised")
        except core.UserFacingError as e:
            good = words.lower() in str(e).lower()
            ok &= good
            print(f"  {'ok  ' if good else 'FAIL'} {label:<34} “{str(e)[:90]}…”")

    print("Construct")
    for n in ("pBAD30.gb", "addgene-plasmid-39296-sequence-49545.gbk", "pFA-KanMX4.dna"):
        check(n, "construct", {"file": f(n)}, lambda r: r["record"]["length"] > 3000)
    rj = load_item(items["pBAD30.gb"])["record"]
    check("digest EcoRI + HindIII", "digest", {"record": rj, "lanes": [["EcoRI", "HindIII"]], "compare_feature": "araC"},
          lambda r: len(r["lanes"][0]["fragments"]) == 2)
    check("enzymes flanking araC", "enzymes", {"record": rj, "filter": "flank", "feature": "araC"}, lambda r: r["upstream"])
    check("suggest digest for araC", "suggest_digest", {"record": rj, "feature": "araC"}, lambda r: r["suggestions"])
    check("ORFs, table 11, alt starts", "orfs", {"record": rj, "table": 11, "min_aa": 60, "starts": "alt"}, lambda r: r["orfs"])
    check("GenBank export", "export", {"record": rj, "format": "genbank"}, lambda r: r["text"].startswith("LOCUS"))
    print("Primers")
    kan = load_item(items["addgene-plasmid-39296-sequence-49545.gbk"])["record"]
    check("KanR pair on pFA6a-kanMX6", "primers",
          {"text": "KanF ATGGGTAAGGAAAAGACTCACG\nKanR TTAGAAAAACTCATCGAGCATC", "preset": "hifi", "template": kan},
          lambda r: r["pcr"]["products"][0]["size"] == 810)
    print("Sanger")
    check("3730.ab1 vs demo reference", "sanger", {"files": [f("3730.ab1")], "reference_file": f("demo_reference.gb"), "demo": True},
          lambda r: r["summary"]["confident"] == 2)
    check("310.ab1 (no quality values)", "sanger", {"files": [f("310.ab1")], "reference": rj},
          lambda r: any("no quality" in x["text"] for x in r["flags"]))
    print("Protein")
    check("KanR from the construct", "protein", {"record": kan, "feature": "KanR"}, lambda r: r["protein"]["ext_reduced"] > 0)
    check("DNA auto-detected and translated", "protein", {"text": "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGATAA"})
    print("Sequence tools")
    check("rotated pBAD30 is the same plasmid", "seqtools",
          {"text": rj["seq"], "circular": True, "other_text": rj["seq"][500:] + rj["seq"][:500], "other_circular": True, "motif": "GAATTC"},
          lambda r: any("Same circular" in x["text"] for x in r["flags"]))
    print("Bad inputs")
    expect_error(".ab1 fed to Construct", "construct", {"file": f("3730.ab1")}, "Sanger trace")
    expect_error("empty file", "construct", {"file": {"name": "empty.gb", "b64": ""}}, "empty")
    expect_error("protein pasted as DNA", "construct", {"text": "MKTAYIAKQRQISFVKSHFSRQLEE"}, "protein")
    expect_error("GenBank fed to Sanger", "sanger", {"files": [f("pBAD30.gb")], "reference": rj}, "not a Sanger trace")
    print("Report and assistant tools")
    try:
        out = tmp / "report.pdf"
        build_pdf(load_item(items["pBAD30.gb"]), out)
        print(f"  ok   PDF report                          {out.stat().st_size // 1024} KB")
    except Exception as e:  # noqa: BLE001
        ok = False
        print(f"  FAIL PDF report: {e}")
    for tool, inp in (("list_items", {}), ("get_item_result", {}), ("digest", {"enzymes": ["EcoRI"]}),
                      ("find_enzymes", {"filter": "outside", "feature": "araC"}), ("translate_feature", {"feature_name": "araC"}),
                      ("find_orfs", {"min_aa": 100}), ("protein_properties", {"seq": "MKTAYIAKQRQISFVKSHFSRQ"}),
                      ("primer_check", {"primers": "a ACGTACGTAGCTAGCTAGG\nb CCTAGCTAGCTACGTACGT"}),
                      ("sanger_summary", {"item": items.get("3730.ab1 vs demo reference")})):
        res = agent.run_tool(tool, inp, items["pBAD30.gb"])
        good = not res.startswith("error") and len(res) > 20
        ok &= good
        print(f"  {'ok  ' if good else 'FAIL'} assistant tool {tool:<20} {len(res):>6} chars")
    print("\nALL PASSED" if ok else "\nSOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
