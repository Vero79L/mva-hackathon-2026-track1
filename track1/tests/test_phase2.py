#!/usr/bin/env python3
"""Phase 2 tests on synthetic data: the pre-registered sensitivity analysis, the aneuploidy scan, and the
report filler. No competition data, no network."""
import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import unittest

from test_track1 import FX1, RULES, SAI1, SCRIPTS, T1, make_work, run

PROTOCOL = os.path.join(T1, "config", "phase2_protocol.yml")
REPORT = os.path.join(T1, "report")


def have(*tools):
    return all(shutil.which(t) for t in tools)


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class Sensitivity(unittest.TestCase):
    def test_runs_every_registered_scenario(self):
        tmp = tempfile.mkdtemp(prefix="t1sens_")
        try:
            w = make_work(os.path.join(tmp, "w"), FX1, SAI1)
            r = run([os.path.join(SCRIPTS, "rank_panel.py"), "--work", w, "--rules", RULES])
            self.assertEqual(r.returncode, 0, r.stderr)
            before = read(os.path.join(w, "trackB_pairs.tsv"))
            out = os.path.join(tmp, "sens.json")
            r = run([os.path.join(SCRIPTS, "sensitivity.py"), "--work", w, "--rules", RULES, "--protocol", PROTOCOL, "--out", out])
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            d = load(out)
            self.assertTrue(d["baseline_reproduced"])
            self.assertEqual(d["n"], 11 * 3 + 4 + 3 + 4 + 4, "11 weights x (half, double, removed) + frequency + SpliceAI + quality + tie-break")
            self.assertEqual(d["n"], d["n_retained"] + d["n_changed"] + d["n_no_pairs"])
            self.assertTrue(all(s["status"] in ("retained", "changed", "no_pairs") for s in d["scenarios"]))
            self.assertEqual(read(os.path.join(w, "trackB_pairs.tsv")), before, "the real outputs must not change")
            self.assertNotIn("15:", json.dumps(d), "no coordinates in the output")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def synth_vcf(path, seed, gain_chr=None, f=0.0, het=2500, hom=1200):
    import numpy as np
    rng = np.random.default_rng(seed)
    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        for c in chroms:
            out.write(f"##contig=<ID={c},length=250000000>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for c in chroms:
            nhet = 0 if c == "Y" else (200 if c == "X" else het)
            pos = np.sort(rng.choice(np.arange(3_000_000, 240_000_000), nhet + hom, replace=False))
            is_het = np.zeros(nhet + hom, bool)
            is_het[rng.choice(nhet + hom, nhet, replace=False)] = True
            scale = 1 + f / 2 if c == gain_chr else 1.0
            for p, h in zip(pos, is_het):
                n = max(1, int(rng.poisson(30 * scale)))
                if h:
                    frac = ((1 + f) / (2 + f) if rng.random() < 0.5 else 1 / (2 + f)) if c == gain_chr else 0.49
                    a, gt = int(rng.binomial(n, frac)), "0/1"
                else:
                    a, gt = n, "1/1"
                out.write(f"{c}\t{p}\t{c}:{p}:A:G\tA\tG\t50\tPASS\t.\tGT:AD\t{gt}:{n - a},{a}\n")
    subprocess.run(["bgzip", "-f", path], check=True)
    return path + ".gz"


@unittest.skipUnless(have("bcftools", "bgzip"), "needs bcftools and bgzip (conda env mva)")
class Aneuploidy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="t1aneu_")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def scan(self, vcf):
        out = os.path.join(self.tmp, os.path.basename(vcf) + ".json")
        r = run([os.path.join(SCRIPTS, "aneuploidy_scan.py"), "--vcf", vcf, "--sex", "MALE", "--protocol", PROTOCOL, "--out", out])
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        return load(out)

    def test_detects_mosaic_gain_of_8(self):
        d = self.scan(synth_vcf(os.path.join(self.tmp, "gain.vcf"), 1, "8", 0.4))
        pre = {x["chrom"]: x["positive"] for x in d["prespecified"]}
        self.assertTrue(pre["8"])
        self.assertFalse(pre["7"])
        self.assertIn({"chrom": "8", "flag": "gain"}, d["exploratory_flags"])

    def test_prespecified_test_finds_a_15_percent_gain_of_7(self):
        d = self.scan(synth_vcf(os.path.join(self.tmp, "gain7.vcf"), 4, "7", 0.15))
        pre = {x["chrom"]: x["positive"] for x in d["prespecified"]}
        self.assertTrue(pre["7"])
        self.assertFalse(pre["8"])

    def test_silent_on_a_normal_genome(self):
        d = self.scan(synth_vcf(os.path.join(self.tmp, "null.vcf"), 2))
        self.assertFalse(any(x["positive"] for x in d["prespecified"]))
        self.assertEqual(d["exploratory_flags"], [])
        self.assertGreater(d["power"]["min_cell_fraction_gain"], 0)
        x = next(r for r in d["chroms"] if r["chrom"] == "X")
        self.assertIsNone(x["z_balance"], "male X is described, not tested")


class Filler(unittest.TestCase):
    def test_fills_template_from_results(self):
        tmp = tempfile.mkdtemp(prefix="t1fill_")
        try:
            w = make_work(os.path.join(tmp, "work"), FX1, SAI1)
            self.assertEqual(run([os.path.join(SCRIPTS, "rank_panel.py"), "--work", w, "--rules", RULES]).returncode, 0)
            res = os.path.join(tmp, "results")
            os.makedirs(os.path.join(res, "report"))
            with open(os.path.join(res, "synth_mva-tiered-v1.csv"), "w", newline="") as fh:
                wr = csv.writer(fh)
                wr.writerow(["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"])
                wr.writerow(["PROBAND01", "chr15", "1000100", "C", "T", "chr15", "1000200", "G", "A", "0.95", "primary", "BUB1B (Tier A); route B #1 pair"])
                wr.writerow(["PROBAND01", "chr3", "3000", "AT", "A", "chr3", "3100", "G", "A", "0.25", "primary", "FANCD2; route A Exomiser exome #1 AR"])
            with open(os.path.join(res, "hpo_terms.tsv"), "w") as fh:
                fh.write("hpo_id\tlabel\tnote\nHP:0000118\tPhenotypic abnormality\t\n")
            with open(os.path.join(res, "report", "sensitivity.json"), "w") as fh:
                json.dump({"n": 2, "n_retained": 2, "scenarios": []}, fh)
            out = os.path.join(res, "report")
            logs = os.path.join(w, "logs")                 # step logs as lib.sh names them; end time = last write
            os.makedirs(logs, exist_ok=True)
            now = time.time()
            stamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(now))
            for name, minutes, text in (("02_fetch_resources", 210, "ok\n"), ("04_qc", 2, "ok\n"),
                                        ("05_annotate", 19, "ok\n"), ("06_exomiser", 60, "\n✗ failed\n")):
                lp = os.path.join(logs, f"{name}_{stamp}.log")
                with open(lp, "w", encoding="utf-8") as fh:
                    fh.write(text)
                os.utime(lp, (now + minutes * 60, now + minutes * 60))
            r = run([os.path.join(REPORT, "fill_report.py"), "--work", w, "--results", res, "--t1", T1, "--out-dir", out,
                     "--no-pdf", "--no-form-in-repo", "--github", "https://github.com/example/x", "--ai-training", "off"])
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            html = read(os.path.join(out, "synth_track1_report.html"))
            m = re.search(r'<script id="report-data" type="application/json">(.*?)</script>', html, re.S)
            self.assertIsNotNone(m)
            d = json.loads(m.group(1).replace("<\\/", "</"))
            self.assertEqual(d["meta"]["team"], "synth")
            self.assertGreater(d["meta"]["n_tests"], 30, "the test count is read from tests/")
            self.assertEqual([a["key"] for a in d["alleles"]], ["15:1000100:C:T", "15:1000200:G:A"])
            self.assertEqual(d["submission"]["rows"][1]["route"], "A")
            self.assertEqual(d["submission"]["rows"][1]["gene"], "FANCD2")
            self.assertEqual(len(d["panel_genes"]), 13)
            self.assertEqual(d["sensitivity"]["n_retained"], 2)
            self.assertIn("version: 1", d["rules_text"])
            self.assertNotIn("SYNTH01", html, "sample names never reach the report")
            self.assertEqual(d["hpo"], {"n": 1}, "only the number of phenotype terms is published")
            self.assertNotIn("HP:0000118", html)
            self.assertFalse({"sex_inferred", "chrx_ratio", "chry_calls"} & set(d["qc"]), "no sex-chromosome metrics")
            self.assertEqual([f["row"] for f in d["form"]], list(range(7, 20)), "all thirteen rows of the organisers' form")
            self.assertTrue(next(f for f in d["form"] if f["row"] == 10)["required"])
            self.assertIn("model training setting off", next(f for f in d["form"] if f["row"] == 10)["a"])
            self.assertLessEqual(d["form_words"], 500)
            self.assertEqual({k: round(v) for k, v in d["runtime"]["steps"].items()}, {"04_qc": 2, "05_annotate": 19},
                             "steps are timed from their logs; a run that stopped with an error is skipped")
            self.assertEqual(round(d["runtime"]["setup_minutes"]), 210)
            for f in d["form"]:
                self.assertNotRegex(f["a"], r"\s,|\(\s*\)|\S  \S", f"row {f['row']}: no gaps left by missing values")
            self.assertTrue(all(p["v1"] in {"15:1000100:C:T", "15:1000200:G:A"} for p in d["routeB"]["pairs_top"]),
                            "only submitted variants are listed")
            for ext in ("md", "xlsx"):
                self.assertTrue(os.path.exists(os.path.join(out, f"synth_track1_methods_form.{ext}")))
            import zipfile
            with zipfile.ZipFile(os.path.join(out, "synth_track1_methods_form.xlsx")) as z:
                sheet = z.read("xl/worksheets/sheet1.xml").decode()
            self.assertIn('r="B10"', sheet)
            self.assertIn("Anthropic Claude", sheet)
            self.assertFalse(os.path.exists(os.path.join(T1, "docs", "synth_track1_methods_form.md")), "tests do not write into the repo")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
