#!/usr/bin/env python3
"""Tests for the generalisation benchmark (make_rules.py, benchmark.py) on synthetic data. No network, no VEP."""
import csv
import gzip
import json
import os
import shutil
import sys
import tempfile
import unittest

import yaml

from test_track1 import RULES, SCRIPTS, T1, run

sys.path.insert(0, SCRIPTS)
import benchmark  # noqa: E402
import make_rules  # noqa: E402

BENCH = os.path.join(SCRIPTS, "benchmark.py")


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def rows(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ann(symbol, cons="intron_variant", impact="MODIFIER", maxaf=".", clnsig=".", clnrev=".", sift=".", polyphen=".",
        mane="NM_1.1", canonical="YES", feature="ENST1"):
    # SYMBOL Feature Consequence IMPACT MANE_SELECT CANONICAL HGVSc HGVSp SIFT PolyPhen gnomADe gnomADg MAX_AF CLNSIG
    # CLNREVSTAT CLNDN Existing
    return [symbol, feature, cons, impact, mane, canonical, ".", ".", sift, polyphen, ".", ".", maxaf, clnsig, clnrev,
            ".", "."]


def site(key, *a, **kw):
    c, p, r, x = key.split(":")
    return [c, p, r, x] + ann(*a, **kw)


GENELISTS = """# synthetic snapshot
panel\tpanelapp_id\tversion\tsymbol\thgnc_id\tconfidence\tmoi
small\t11\t1.0\tGENEA\tHGNC:1\t3\tBIALLELIC, autosomal or pseudoautosomal
small\t11\t1.0\tGENEB\tHGNC:2\t3\tMONOALLELIC, autosomal or pseudoautosomal, NOT imprinted
small\t11\t1.0\tGENEC\tHGNC:3\t2\tBIALLELIC, autosomal or pseudoautosomal
small\t11\t1.0\tGENER\tHGNC:4\t1\tBIALLELIC, autosomal or pseudoautosomal
small\t11\t1.0\tGENEX\tHGNC:5\t3\tX-LINKED: hemizygous mutation in males, biallelic mutations in females
small\t11\t1.0\tNOGTF\tHGNC:6\t3\tBIALLELIC, autosomal or pseudoautosomal
big\t22\t2.0\tGENEA\tHGNC:1\t3\tBIALLELIC, autosomal or pseudoautosomal
big\t22\t2.0\tGENED\tHGNC:7\t3\tBOTH monoallelic and biallelic, autosomal or pseudoautosomal
"""

GTF = [("15", 1000000, 1010000, "BUB1B"), ("11", 2000000, 2005000, "CEP57"), ("2", 3000, 9000, "BUB1"),
       ("7", 100000, 110000, "GENEA"), ("7", 200000, 205000, "GENEB"), ("8", 50000, 60000, "GENEC"),
       ("8", 70000, 80000, "GENER"), ("X", 5000, 9000, "GENEX"), ("9", 400000, 420000, "GENED")]


def protocol(path, samples, n=3):
    with open(os.path.join(T1, "config", "phase2_protocol.yml"), encoding="utf-8") as fh:
        p = yaml.safe_load(fh)
    p["benchmark"]["backgrounds"]["samples"] = samples
    p["benchmark"]["panels"]["panelapp"] = {"small": {"id": 11, "version": "1.0"}, "big": {"id": 22, "version": "2.0"}}
    p["benchmark"]["spikes"]["per_panel_per_background"] = n
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(p, fh, sort_keys=False, allow_unicode=True)
    return path


def write_gtf(path):
    with gzip.open(path, "wt") as fh:
        for c, s, e, g in GTF:
            fh.write(f'{c}\tensembl\tgene\t{s}\t{e}\t.\t+\t.\tgene_id "ENSG{g}"; gene_name "{g}"; gene_biotype "protein_coding";\n')
            fh.write(f'{c}\tensembl\texon\t{s}\t{s + 100}\t.\t+\t.\tgene_id "ENSG{g}"; gene_name "{g}";\n')


class Rules(unittest.TestCase):
    def test_swap_changes_only_the_knowledge_layer(self):
        tmp = tempfile.mkdtemp(prefix="t1mr_")
        try:
            gl = os.path.join(tmp, "g.tsv")
            with open(gl, "w") as fh:
                fh.write(GENELISTS)
            lists = make_rules.read_genelists(gl)
            self.assertEqual(lists["small"]["tier_a"], ["GENEA", "GENEB", "GENEX", "NOGTF"])
            self.assertEqual(lists["small"]["tier_b"], ["GENEC"], "amber is tier B, red is not used")
            self.assertIn("GENED", lists["big"]["biallelic"], "'BOTH monoallelic and biallelic' counts as recessive")
            self.assertNotIn("GENEB", lists["small"]["biallelic"])
            out = os.path.join(tmp, "r.yml")
            r = run([os.path.join(SCRIPTS, "make_rules.py"), "--genelists", gl, "--panel", "small", "--out", out])
            self.assertEqual(r.returncode, 0, r.stderr)
            base = make_rules.load_rules(RULES)
            new = make_rules.load_rules(out)
            self.assertEqual(make_rules.engine_part(new), make_rules.engine_part(base), "weights and thresholds unchanged")
            self.assertEqual(new["genes"]["tier_a"], lists["small"]["tier_a"])
            self.assertEqual(new["known_alleles"], [])
            r = run([os.path.join(SCRIPTS, "make_rules.py"), "--tier-a", "X1,X2", "--tier-b", "X2,X3", "--out", out])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(make_rules.load_rules(out)["genes"]["tier_b"], ["X3"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_sites_format_matches_the_pipeline_columns(self):
        fmt = benchmark.sites_format(os.path.join(T1, "config", "splitvep_format.txt"))
        toks = fmt.replace("\\n", "").split("\\t")
        self.assertEqual(toks[:4], ["%CHROM", "%POS", "%REF", "%ALT"])
        self.assertEqual(len(toks), 4 + len(benchmark.ANN))
        self.assertEqual(toks[4], "%SYMBOL")
        self.assertEqual(toks[-1], "%Existing_variation")


class Pipeline(unittest.TestCase):
    """regions -> spikes -> catalog -> simulate on a synthetic world with the real rank_panel.py."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="t1bench_")
        t = cls.tmp
        cls.proto = protocol(os.path.join(t, "proto.yml"), ["S1", "S2"])
        cls.gl = os.path.join(t, "genelists.tsv")
        with open(cls.gl, "w") as fh:
            fh.write(GENELISTS)
        write_gtf(os.path.join(t, "g.gtf.gz"))
        cls.reg = os.path.join(t, "regions")
        cls.r_regions = run([BENCH, "regions", "--gtf", os.path.join(t, "g.gtf.gz"), "--genelists", cls.gl,
                             "--rules", RULES, "--protocol", cls.proto, "--flank", "1000", "--out", cls.reg])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_1_regions(self):
        self.assertEqual(self.r_regions.returncode, 0, self.r_regions.stderr)
        P = json.loads(read(os.path.join(self.reg, "panels.json")))["panels"]
        self.assertEqual(P["small"]["spike_genes"], ["GENEA"], "green, recessive and autosomal only")
        self.assertEqual(P["small"]["missing_from_gtf"], ["NOGTF"])
        self.assertEqual(P["big"]["spike_genes"], ["GENEA", "GENED"])
        self.assertEqual(P["mva"]["rules"], os.path.abspath(RULES), "MVA uses the frozen file itself")
        self.assertIn("15\t40117037\t40117138\tBUB1B:known", read(os.path.join(self.reg, "mva.bed")))
        union = read(os.path.join(self.reg, "union.bed"))
        self.assertNotIn("X\t", union, "backgrounds are autosomal, so the shared regions are too")
        self.assertTrue(read(os.path.join(self.reg, "union.chr.bed")).startswith("chr"))
        before = os.path.getmtime(os.path.join(self.reg, "union.bed"))
        r = run([BENCH, "regions", "--gtf", os.path.join(self.tmp, "g.gtf.gz"), "--genelists", self.gl, "--rules", RULES,
                 "--protocol", self.proto, "--flank", "1000", "--out", self.reg])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(os.path.getmtime(os.path.join(self.reg, "union.bed")), before, "unchanged outputs are not rewritten")

    def test_2_spikes(self):
        raw = os.path.join(self.tmp, "raw.tsv")
        lines = [
            "7\t100100\tC\tT\tPathogenic\tcriteria_provided,_multiple_submitters,_no_conflicts\tGENEA:1\t1",
            "7\t100200\tG\tA\tLikely_pathogenic\tcriteria_provided,_single_submitter\tGENEA:1\t1",
            "7\t100300\tG\tA\tPathogenic\tno_assertion_criteria_provided\tGENEA:1\t1",          # 0 stars
            "7\t100400\tG\tA\tConflicting_classifications_of_pathogenicity\tcriteria_provided,_conflicting_classifications\tGENEA:1\t1",
            "7\t100500\tG\tA\tPathogenic|risk_factor\tcriteria_provided,_single_submitter\tGENEA:1\t2",  # somatic only
            "7\t200100\tC\tT\tPathogenic\tcriteria_provided,_single_submitter\tGENEB:2\t1",          # dominant gene
            "9\t400100\tC\tT\tPathogenic/Likely_pathogenic\treviewed_by_expert_panel\tGENED:7\t1",
            "9\t400200\t" + "A" * 60 + "\tA\tPathogenic\tcriteria_provided,_single_submitter\tGENED:7\t1",  # too long
            "X\t6000\tC\tT\tPathogenic\tcriteria_provided,_single_submitter\tGENEX:5\t1",            # not autosomal
        ]
        with open(raw, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        out = os.path.join(self.tmp, "spikes.tsv")
        r = run([BENCH, "spikes", "--raw", raw, "--panels", os.path.join(self.reg, "panels.json"), "--out", out,
                 "--vcf", os.path.join(self.tmp, "spikes.vcf")])
        self.assertEqual(r.returncode, 0, r.stderr)
        got = {x["key"]: x for x in rows(out)}
        self.assertEqual(sorted(got), ["7:100100:C:T", "7:100200:G:A", "9:400100:C:T"])
        self.assertEqual(got["7:100200:G:A"]["stars"], "1")
        self.assertIn("7\t100100\t7:100100:C:T\tC\tT", read(os.path.join(self.tmp, "spikes.vcf")))

    def test_3_catalog_and_simulation(self):
        t = self.tmp
        A1, A2, A3 = "7:100100:C:T", "7:100200:G:A", "7:100250:G:C"         # GENEA pathogenic (A3 common)
        D1, D2 = "9:400100:C:T", "9:400150:A:G"
        BG_LOF1, BG_LOF2 = "9:401000:C:T", "9:402000:G:A"                   # background LoF pair in GENED
        COMMON = "7:105000:A:G"
        sites = [
            site(A1, "GENEA", "stop_gained", "HIGH", clnsig="Pathogenic", clnrev="criteria_provided&_multiple_submitters&_no_conflicts"),
            site(A1, "GENEA", "intron_variant", mane=".", canonical=".", feature="ENST2", clnsig="Pathogenic"),
            site(A2, "GENEA", "missense_variant", "MODERATE", sift="deleterious(0)", polyphen="probably_damaging(1)",
                 clnsig="Likely_pathogenic", clnrev="criteria_provided&_single_submitter"),
            site(A3, "GENEA", "missense_variant", "MODERATE", maxaf="0.2", clnsig="Pathogenic",
                 clnrev="criteria_provided&_single_submitter"),
            site(D1, "GENED", "missense_variant", "MODERATE", sift="deleterious(0)", polyphen="probably_damaging(1)",
                 clnsig="Pathogenic", clnrev="reviewed_by_expert_panel"),
            site(D2, "GENED", "missense_variant", "MODERATE", sift="deleterious(0)", polyphen="probably_damaging(1)",
                 clnsig="Pathogenic", clnrev="reviewed_by_expert_panel"),
            site(BG_LOF1, "GENED", "stop_gained", "HIGH", maxaf="0.0001"),
            site(BG_LOF2, "GENED", "frameshift_variant", "HIGH", maxaf="0.0002"),
            site(COMMON, "GENEA", "intron_variant", maxaf="0.3"),
            site("15:40117088:G:A", "BMF", "upstream_gene_variant", maxaf="0.02"),
            site("1:500:A:G", "OTHERGENE", "stop_gained", "HIGH"),
        ]
        raw = os.path.join(t, "splitvep.tsv")
        with open(raw, "w") as fh:
            fh.write("".join("\t".join(s) + "\n" for s in sites))
        cat = os.path.join(t, "catalog.tsv.gz")
        r = run([BENCH, "catalog", "--panels", os.path.join(self.reg, "panels.json"), "--rules", RULES, "--input", raw,
                 "--out", cat])
        self.assertEqual(r.returncode, 0, r.stderr)
        with gzip.open(cat, "rt") as fh:
            kept = {benchmark.variant_key(*ln.split("\t")[:4]) for ln in fh}
        self.assertNotIn(COMMON, kept, "common and not pathogenic: can never qualify")
        self.assertNotIn("1:500:A:G", kept, "not in any panel")
        self.assertIn("15:40117088:G:A", kept, "known allele kept although intergenic and 2% frequent")
        self.assertIn(A3, kept, "ClinVar pathogenic kept even when common")

        spikes = os.path.join(t, "spikes2.tsv")
        with open(spikes, "w") as fh:
            fh.write("key\tgenes\tclnsig\tstars\n")
            for k, g in ((A1, "GENEA"), (A2, "GENEA"), (A3, "GENEA"), (D1, "GENED"), (D2, "GENED")):
                fh.write(f"{k}\t{g}\tPathogenic\t2\n")
        gdir = os.path.join(t, "bg")
        os.makedirs(gdir, exist_ok=True)
        for s, extra in (("S1", [BG_LOF1, BG_LOF2]), ("S2", [])):
            with open(os.path.join(gdir, s + ".gt.tsv"), "w") as fh:
                for k in [COMMON] + extra:
                    c, p, ref, alt = k.split(":")
                    fh.write(f"{c}\t{p}\t{ref}\t{alt}\t{k}\tPASS\t50\t0/1\t20,20\t40\t99\n")
        out = os.path.join(t, "results")
        r = run([BENCH, "simulate", "--regions", self.reg, "--catalog", cat, "--genotypes", gdir, "--spikes", spikes,
                 "--rules", RULES, "--protocol", self.proto, "--out", out, "--workers", "2", "--tmp", os.path.join(t, "tmp"),
                 "--protocol-commit", "abc1234"])
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        d = json.loads(read(os.path.join(out, "benchmark.json")))
        P = {p["key"]: p for p in d["panels"]}
        self.assertEqual(set(P), {"mva", "small", "big"})
        self.assertEqual(P["small"]["n_eligible_variants"], 2, "the common pathogenic variant is not a spike")
        self.assertEqual(P["small"]["modes"]["clinvar"]["rank1"], 1.0)
        self.assertEqual(P["small"]["modes"]["clinvar"]["n"], 2 * 3)
        self.assertEqual(P["mva"]["modes"]["novel"]["n"], 0, "no eligible MVA spikes in this synthetic world")
        self.assertEqual(P["big"]["healthy"]["per_background"]["S1"]["candidates"], 1, "S1 carries a background LoF pair")
        self.assertEqual(P["big"]["healthy"]["per_background"]["S2"]["candidates"], 0)
        sims = rows(os.path.join(out, "simulations.tsv"))
        # S1 carries a rare LoF pair in GENED (9 points each). Expected ranks under the frozen rules, novel mode:
        #   GENEA LoF/missense: pair scores min(10, 4) = 4, behind the background pair -> 2
        #   GENED missense/missense: 4 points, tied with the four mixed (background LoF, spiked missense) pairs,
        #     which win the sum_points tie-break (13 > 8) -> 6
        #   homozygous LoF: 10 points, ahead of the background pair -> 1; homozygous missense: 4 points -> 2
        expect = {"LoF/missense": "2", "missense/missense": "6", "hom LoF": "1", "hom missense": "2"}
        s1 = [x for x in sims if x["panel"] == "big" and x["background"] == "S1" and x["mode"] == "novel"]
        self.assertTrue(s1)
        for x in s1:
            self.assertEqual(x["rank"], expect[x["cls"]], x)
        s2 = [x for x in sims if x["panel"] in ("big", "small") and x["background"] == "S2" and x["mode"] != "healthy"]
        self.assertTrue(s2 and all(x["rank"] == "1" for x in s2), "with no competing candidate every genotype ranks first")
        self.assertEqual(d["protocol_commit"], "abc1234")
        self.assertNotIn("PROBAND", json.dumps(d))


if __name__ == "__main__":
    unittest.main(verbosity=2)
