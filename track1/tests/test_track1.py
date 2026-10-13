#!/usr/bin/env python3
"""Track 1 的測試：全部用合成資料（虛構的變異與表現型），不需要比賽資料、bcftools、VEP 或 Exomiser。

執行：bash track1/tests/run_tests.sh
涵蓋：排序規則（v1.2 的決勝順序、cis 排除、1/2 拆分、已知低表現等位基因、SpliceAI 基因比對、舊基因名稱）、
提交 CSV（merge 順序、EPCR、chr 前綴、Exomiser 寫法換回 VCF 寫法、格式錯誤偵測、評分模擬、官方介面轉接）、
HPO 整理、phenopacket、panel 範圍、共用 shell 函式（含下載續傳與換版）、事先登錄檢查、
macOS 相容性（含 $var 後面直接接中文）。
"""
import csv
import gzip
import http.server
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
T1 = os.path.dirname(HERE)
SCRIPTS = os.path.join(T1, "scripts")
RULES = os.path.join(T1, "config", "ranking_rules.yml")
sys.path.insert(0, SCRIPTS)
import t1common  # noqa: E402

PY = sys.executable
COLS = 30


def key(v):
    return "%s:%d:%s:%s" % v


def vep_row(v, gt, symbol, cons, impact, *, vid=None, filt="PASS", gq="99", dp="35", pgt=".", pid=".",
            mane="NM_000001.1", canonical="YES", hgvsc=".", hgvsp=".", sift=".", polyphen=".",
            maxaf=".", clnsig=".", clnrev=".", feature="ENST00000000001"):
    chrom, pos, ref, alt = v
    r = [chrom, str(pos), vid or key(v), ref, alt, filt, "50", gt, "15,15", dp, gq, pgt, pid,
         symbol, feature, cons, impact, mane, canonical, hgvsc, hgvsp, sift, polyphen, ".", ".", maxaf,
         clnsig, clnrev, ".", "."]
    assert len(r) == COLS
    return "\t".join(r)


def write_vcf(path, variants, sample="SYNTH01"):
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=2,length=242193529>", "##contig=<ID=3,length=198295559>",
             "##contig=<ID=9,length=138394717>", "##contig=<ID=11,length=135086622>",
             "##contig=<ID=15,length=101991189>", "##contig=<ID=17,length=83257441>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample]
    order = {"2": 2, "3": 3, "9": 9, "11": 11, "15": 15, "17": 17}
    for c, p, r, a in sorted(set(variants), key=lambda x: (order[x[0]], x[1], x[3])):
        lines.append(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\tGT\t0/1")
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def run(args, **kw):
    return subprocess.run([PY] + args, capture_output=True, text=True, **kw)


# 合成情境 1：BUB1B 有截斷（T）、罕見錯義（M）、和 T 證實 trans 的罕見內含子（I）、常見內含子（C）、
# SpliceAI 中度的同義（S）、非 PASS 錯義（N）、和 T 證實 cis 的內含子（Q）；CEP57 只有一個截斷；BUB1 是 1/2 拆分。
T = ("15", 1000100, "C", "T")
M = ("15", 1000200, "G", "A")
I = ("15", 1000150, "A", "G")  # noqa: E741
C = ("15", 1000300, "T", "C")
S = ("15", 1000400, "C", "G")
N = ("15", 1000600, "G", "C")
Q = ("15", 1000160, "G", "T")
CEP = ("11", 2000100, "C", "A")
B1a, B1b = ("2", 5000, "A", "C"), ("2", 5000, "A", "G")
OTHER = ("17", 7000, "C", "T")
GX1, GX2 = ("3", 3000, "AT", "A"), ("3", 3100, "G", "A")
SF = ("17", 9000, "G", "A")
PID = "1000100_C_T"

FX1 = [
    vep_row(T, "0|1", "BUB1B", "stop_gained", "HIGH", pgt="0|1", pid=PID, hgvsp="ENSP1:p.Arg100Ter"),
    vep_row(T, "0|1", "BUB1B", "intron_variant", "MODIFIER", mane=".", canonical=".", feature="ENST00000000002"),
    vep_row(M, "0/1", "BUB1B", "missense_variant", "MODERATE", sift="tolerated(0.21)", polyphen="probably_damaging(0.97)",
            hgvsp="ENSP1:p.Gly200Arg"),
    vep_row(I, "0|1", "BUB1B", "intron_variant", "MODIFIER", pgt="1|0", pid=PID, maxaf="5e-05"),
    vep_row(C, "0/1", "BUB1B", "intron_variant", "MODIFIER", maxaf="0.21"),
    vep_row(S, "0/1", "BUB1B", "synonymous_variant", "LOW", maxaf="0.003"),
    vep_row(N, "0/1", "BUB1B", "missense_variant", "MODERATE", filt="LowQual", gq="12",
            sift="deleterious(0.01)", polyphen="benign(0.1)"),
    vep_row(Q, "0|1", "BUB1B", "intron_variant", "MODIFIER", pgt="0|1", pid=PID),
    vep_row(CEP, "0/1", "CEP57", "stop_gained", "HIGH"),
    vep_row(B1a, "1/0", "BUB1", "missense_variant", "MODERATE", vid="2:5000:A:C,G", sift="deleterious(0)",
            polyphen="probably_damaging(1)"),
    vep_row(B1b, "0/1", "BUB1", "missense_variant", "MODERATE", vid="2:5000:A:C,G", sift="deleterious(0)",
            polyphen="probably_damaging(1)"),
    vep_row(OTHER, "0/1", "TP53", "stop_gained", "HIGH"),
]
SAI1 = [f"15\t{S[1]}\tC\tG\tG|BUB1B|0.35|0.00|0.02|0.00|-3|10|5|-20",
        f"15\t{M[1]}\tG\tA\tA|BUB1B|0.01|0.00|0.00|0.00|1|2|3|4"]

# 合成情境 2：BUB1B 截斷（T）+ ClinGen 已知的上游低表現等位基因（H，基因間、頻率 2%）+ 有害錯義（M2）；
# H 的 SpliceAI 註解屬於別的基因（BMF），不能算進 BUB1B；SMC5 有 ClinVar 1 星；KNL1 用舊名 CASC5。
H = ("15", 40117088, "G", "A")
M2 = ("15", 1000200, "G", "A")
SMC = ("9", 3000100, "C", "G")
KN = ("15", 2000500, "A", "T")
FX2 = [
    vep_row(T, "0/1", "BUB1B", "stop_gained", "HIGH"),
    vep_row(H, "0/1", "BMF", "upstream_gene_variant", "MODIFIER", maxaf="0.02"),
    vep_row(M2, "0/1", "BUB1B", "missense_variant", "MODERATE", sift="deleterious(0.01)", polyphen="probably_damaging(0.99)"),
    vep_row(SMC, "0/1", "SMC5", "missense_variant", "MODERATE", maxaf="0.0005", clnsig="Likely_pathogenic",
            clnrev="criteria_provided&_single_submitter"),
    vep_row(KN, "0/1", "CASC5", "missense_variant", "MODERATE", sift="deleterious(0)", polyphen="possibly_damaging(0.8)"),
]
SAI2 = [f"15\t{H[1]}\tG\tA\tA|BMF|0.90|0.00|0.00|0.00|1|1|1|1"]


def make_work(tmp, rows, sai=None):
    os.makedirs(tmp, exist_ok=True)
    with open(os.path.join(tmp, "panel.vep.tsv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    if sai:
        with open(os.path.join(tmp, "panel.spliceai.tsv"), "w") as fh:
            fh.write("\n".join(sai) + "\n")
    return tmp


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class RankPanel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="t1test_")
        cls.w1 = make_work(os.path.join(cls.tmp, "fx1"), FX1, SAI1)
        cls.r1 = run([os.path.join(SCRIPTS, "rank_panel.py"), "--work", cls.w1, "--rules", RULES])
        cls.w2 = make_work(os.path.join(cls.tmp, "fx2"), FX2, SAI2)
        cls.r2 = run([os.path.join(SCRIPTS, "rank_panel.py"), "--work", cls.w2, "--rules", RULES])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_runs(self):
        self.assertEqual(self.r1.returncode, 0, self.r1.stderr + self.r1.stdout)
        self.assertEqual(self.r2.returncode, 0, self.r2.stderr + self.r2.stdout)

    def pairs(self, w):
        return read_tsv(os.path.join(w, "trackB_pairs.tsv"))

    def variants(self, w):
        return {(r["key"], r["gene"]): r for r in read_tsv(os.path.join(w, "trackB_variants.tsv"))}

    def test_rank1_missense_beats_phased_intronic(self):
        p = self.pairs(self.w1)
        self.assertEqual({p[0]["v1"], p[0]["v2"]}, {key(T), key(M)}, "第 1 名應是截斷 + 罕見錯義")
        self.assertEqual({p[1]["v1"], p[1]["v2"]}, {key(T), key(N)}, "同分時 PASS 的配對排前面")
        self.assertEqual({p[3]["v1"], p[3]["v2"]}, {key(T), key(I)}, "trans 的內含子配對排在錯義配對之後")

    def test_mane_transcript_used(self):
        v = self.variants(self.w1)
        self.assertIn("PVS1", v[(key(T), "BUB1B")]["evidence"])
        self.assertEqual(v[(key(T), "BUB1B")]["points"], "10")

    def test_common_not_qualifying(self):
        self.assertEqual(self.variants(self.w1)[(key(C), "BUB1B")]["qualifying"], "False")

    def test_spliceai_supporting(self):
        r = self.variants(self.w1)[(key(S), "BUB1B")]
        self.assertIn("SpliceAI_supporting", r["evidence"])
        self.assertEqual(r["qualifying"], "True")

    def test_cis_excluded(self):
        ex = read_tsv(os.path.join(self.w1, "trackB_excluded.tsv"))
        self.assertTrue(any({e["v1"], e["v2"]} == {key(T), key(Q)} for e in ex))
        self.assertFalse(any({p["v1"], p["v2"]} == {key(T), key(Q)} for p in self.pairs(self.w1)))
        iq = [p for p in self.pairs(self.w1) if {p["v1"], p["v2"]} == {key(I), key(Q)}]
        self.assertEqual(iq[0]["phase"], "trans")

    def test_split_multiallelic_is_trans(self):
        b = [p for p in self.pairs(self.w1) if p["gene"] == "BUB1"]
        self.assertEqual(len(b), 1)
        self.assertEqual(b[0]["phase"], "trans")
        self.assertEqual(b[0]["tier"], "B")

    def test_single_het_after_pairs(self):
        p = self.pairs(self.w1)
        cep = [x for x in p if x["gene"] == "CEP57"][0]
        self.assertEqual(cep["kind"], "single_het_second_allele_missing")
        self.assertEqual(float(cep["pair_points"]), 5.0)
        last_a_pair = max(int(x["rank"]) for x in p if x["tier"] == "A" and x["kind"] == "pair")
        self.assertGreater(int(cep["rank"]), last_a_pair)
        self.assertLess(int(cep["rank"]), int([x for x in p if x["gene"] == "BUB1"][0]["rank"]))

    def test_non_panel_gene_ignored(self):
        self.assertFalse(any(k[1] == "TP53" for k in self.variants(self.w1)))

    def test_known_hypomorph_rank1(self):
        p = self.pairs(self.w2)
        self.assertEqual({p[0]["v1"], p[0]["v2"]}, {key(T), key(H)}, "截斷 + 已知低表現等位基因應排第 1")
        h = self.variants(self.w2)[(key(H), "BUB1B")]
        self.assertEqual(h["qualifying"], "True", "已知等位基因不受 1% 頻率門檻限制")
        self.assertIn("Known_allele", h["evidence"])
        self.assertNotIn("SpliceAI", h["evidence"], "別的基因（BMF）的 SpliceAI 分數不能算進 BUB1B")

    def test_clinvar_and_alias(self):
        v = self.variants(self.w2)
        self.assertIn("ClinVar_PLP_1star", v[(key(SMC), "SMC5")]["evidence"])
        self.assertIn((key(KN), "KNL1"), v, "舊名 CASC5 應換成 KNL1")


class Submission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="t1sub_")
        cls.work = make_work(os.path.join(cls.tmp, "work"), FX1, SAI1)
        cls.res = os.path.join(cls.tmp, "results")
        os.makedirs(cls.res)
        r = run([os.path.join(SCRIPTS, "rank_panel.py"), "--work", cls.work, "--rules", RULES])
        assert r.returncode == 0, r.stderr
        allv = [T, M, I, C, S, N, Q, CEP, B1a, B1b, GX1, GX2, SF]
        write_vcf(os.path.join(cls.work, "norm.vcf.gz"), allv)
        ex = os.path.join(cls.work, "exomiser")
        os.makedirs(ex)
        hdr = ["#RANK", "ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "MOI", "P-VALUE", "EXOMISER_GENE_COMBINED_SCORE",
               "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_VARIANT_SCORE", "EXOMISER_VARIANT_SCORE",
               "CONTRIBUTING_VARIANT", "WHITELIST_VARIANT", "VCF_ID", "RS_ID", "CONTIG", "START", "END", "REF", "ALT",
               "CHANGE_LENGTH", "QUAL", "FILTER", "GENOTYPE"]
        rows = [
            ["1", "x", "GENEX", "1", "AR", "0", "0.9", "0.8", "0.9", "0.9", "1", "0", ".", ".", "3", "3001", "3001", "T", "", "-1", "50", "PASS", "0/1"],
            ["1", "x", "GENEX", "1", "AR", "0", "0.9", "0.8", "0.9", "0.9", "1", "0", ".", ".", "3", "3100", "3100", "G", "A", "0", "50", "PASS", "0/1"],
            ["2", "x", "GENEY", "2", "AD", "0", "0.8", "0.8", "0.8", "0.8", "1", "0", ".", ".", "17", "9000", "9000", "G", "A", "0", "50", "PASS", "0/1"],
            ["3", "x", "BUB1B", "3", "AR", "0", "0.7", "0.7", "0.7", "0.7", "1", "0", ".", ".", "15", str(T[1]), str(T[1]), "C", "T", "0", "50", "PASS", "0|1"],
            ["3", "x", "BUB1B", "3", "AR", "0", "0.7", "0.7", "0.7", "0.7", "1", "0", ".", ".", "15", str(M[1]), str(M[1]), "G", "A", "0", "50", "PASS", "0/1"],
        ]
        with open(os.path.join(ex, "PROBAND01_exome.variants.tsv"), "w") as fh:
            fh.write("\t".join(hdr) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
        with open(os.path.join(cls.res, "secondary_findings.tsv"), "w") as fh:
            fh.write("chrom\tpos\tref\talt\tgene\treason\nchr17\t9000\tG\tA\tGENEZ\tsynthetic test\n")
        cls.draft = run([os.path.join(SCRIPTS, "make_submission.py"), "draft", "--work", cls.work, "--results", cls.res,
                         "--rules", RULES, "--name", "tester"])
        cls.csv = os.path.join(cls.res, "tester_mva-tiered-v1.csv")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def rows(self):
        with open(self.csv, newline="") as fh:
            return list(csv.DictReader(fh))

    def check(self, path, scorer=None):
        args = [os.path.join(SCRIPTS, "make_submission.py"), "check", "--csv", path,
                "--vcf", os.path.join(self.work, "norm.vcf.gz"), "--scorer", scorer or "/nonexistent"]
        return run(args)

    def test_draft_layout(self):
        self.assertEqual(self.draft.returncode, 0, self.draft.stderr + self.draft.stdout)
        rows = self.rows()
        self.assertEqual(len(rows), 9)
        r1 = rows[0]
        self.assertEqual({(r1["chrom_1"], r1["pos_1"]), (r1["chrom_2"], r1["pos_2"])},
                         {("chr15", str(T[1])), ("chr15", str(M[1]))})
        self.assertEqual(r1["epcr"], "0.95")
        self.assertTrue(all("BUB1B" in r["notes"] for r in rows[1:4]), "第 2–4 列是同一基因的其他配對")
        self.assertIn("GENEX", rows[4]["notes"], "第 5 列來自 Exomiser 的隱性基因")
        self.assertEqual((rows[4]["chrom_1"], rows[4]["pos_1"], rows[4]["ref_1"], rows[4]["alt_1"]),
                         ("chr3", "3000", "AT", "A"), "Exomiser 的 indel 寫法要換回標準化 VCF 的寫法")
        self.assertFalse(any("GENEY" in r["notes"] for r in rows), "顯性模式的基因不放進配對")
        self.assertEqual(rows[-1]["finding_type"], "secondary")
        self.assertEqual(rows[-1]["epcr"], "0.05")
        eps = [float(r["epcr"]) for r in rows]
        self.assertEqual(len(set(eps)), len(eps))
        self.assertEqual(eps, sorted(eps, reverse=True))
        sets = [frozenset(((r["chrom_1"], r["pos_1"]), (r["chrom_2"], r["pos_2"]))) for r in rows]
        self.assertEqual(len(set(sets)), len(sets), "不能有重複的列")

    def test_check_passes_and_simulates(self):
        r = self.check(self.csv)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("格式檢查通過", r.stdout)
        self.assertIn("答案 = 第 1 名（EPCR 0.95）→ Rank points 100、F-max 1.000", r.stdout)
        self.assertIn("第 1 名只對一個變異（另一個不在清單裡）→ Rank points 50", r.stdout)

    def test_official_interface_adapter(self):
        r = self.check(self.csv, scorer=os.path.join(HERE, "mock_evaluation.py"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("官方 Rank points 100", r.stdout)
        self.assertIn("和本機一致", r.stdout)

    def mutate(self, fn, name):
        rows = self.rows()
        fn(rows)
        p = os.path.join(self.tmp, name)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        return self.check(p)

    def test_detects_missing_chr(self):
        def f(rows):
            rows[0]["chrom_1"] = "15"
        r = self.mutate(f, "nochr.csv")
        self.assertEqual(r.returncode, 1)
        self.assertIn("chr15", r.stdout)

    def test_detects_duplicate_epcr(self):
        def f(rows):
            rows[1]["epcr"] = rows[0]["epcr"]
        r = self.mutate(f, "dup.csv")
        self.assertEqual(r.returncode, 1)
        self.assertIn("EPCR 有重複", r.stdout)

    def test_detects_wrong_allele(self):
        def f(rows):
            rows[0]["alt_1"] = "G" if rows[0]["alt_1"] != "G" else "C"
        r = self.mutate(f, "wrongalt.csv")
        self.assertEqual(r.returncode, 1)
        self.assertIn("找不到", r.stdout)

    def test_detects_too_many_rows(self):
        def f(rows):
            for k in range(3):
                extra = dict(rows[-1])
                extra["epcr"] = "0.0%d" % (3 - k)
                rows.append(extra)
        r = self.mutate(f, "toomany.csv")
        self.assertEqual(r.returncode, 1)
        self.assertIn("1–10 列", r.stdout)

    def test_simulate_semantics(self):
        sys.path.insert(0, SCRIPTS)
        import make_submission as ms
        a, b, c = ("chr1", 1, "A", "C"), ("chr1", 2, "G", "T"), ("chr1", 3, "T", "A")
        rows = [dict(i=1, variants=frozenset({a, b}), epcr=0.9), dict(i=2, variants=frozenset({a, c}), epcr=0.5)]
        self.assertEqual(ms.simulate(rows, frozenset({a, b}))[:2], (100, 1.0))
        self.assertEqual(ms.simulate(rows, frozenset({a, c}))[0], 50)
        tie = [dict(i=1, variants=frozenset({a, b}), epcr=0.9), dict(i=2, variants=frozenset({c}), epcr=0.9)]
        pts, fm, _, _ = ms.simulate(tie, frozenset({a, b}))
        self.assertEqual(pts, 100)
        self.assertLess(fm, 1.0, "和第 1 名同分的列會拉低 F-max")


class HpoAndPhenopacket(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="t1hpo_")
        self.obo = os.path.join(self.tmp, "hp.obo")
        with open(self.obo, "w") as fh:
            fh.write("format-version: 1.2\ndata-version: hp/releases/2099-01-01\n\n"
                     "[Term]\nid: HP:0000001\nname: All\n\n"
                     "[Term]\nid: HP:0000118\nname: Phenotypic abnormality\nalt_id: HP:0000005\n\n"
                     "[Term]\nid: HP:9999999\nname: obsolete synthetic term\nis_obsolete: true\nreplaced_by: HP:0000118\n")
        self.data = os.path.join(self.tmp, "data")
        os.makedirs(self.data)
        with open(os.path.join(self.data, "phenotype.txt"), "w") as fh:
            fh.write("Synthetic notes\nHP:0000001 All\nHP:9999999 - old term\nHP:0000005: alt id\n")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_extract_and_phenopacket(self):
        out = os.path.join(self.tmp, "hpo_terms.tsv")
        r = run([os.path.join(SCRIPTS, "hpo_terms.py"), "extract", "--data", self.data, "--out", out, "--obo", self.obo])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        ids = [x["hpo_id"] for x in read_tsv(out)]
        self.assertEqual(ids, ["HP:0000001", "HP:0000118"], "過時與舊編號都換成 HP:0000118，且不重複")
        vcf = os.path.join(self.tmp, "x.vcf.gz")
        write_vcf(vcf, [T])
        pp = os.path.join(self.tmp, "pp.yml")
        r = run([os.path.join(SCRIPTS, "make_phenopacket.py"), "--hpo", out, "--vcf", vcf, "--sex", "female",
                 "--obo", self.obo, "--out", pp])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        import yaml
        d = yaml.safe_load(open(pp))
        self.assertEqual(d["subject"], {"id": "SYNTH01", "sex": "FEMALE"})
        self.assertEqual([f["type"]["id"] for f in d["phenotypicFeatures"]], ["HP:0000001", "HP:0000118"])
        self.assertEqual(d["metaData"]["phenopacketSchemaVersion"], "2.0")

    def test_no_hpo_found(self):
        empty = os.path.join(self.tmp, "empty")
        os.makedirs(empty)
        r = run([os.path.join(SCRIPTS, "hpo_terms.py"), "extract", "--data", empty, "--out",
                 os.path.join(self.tmp, "h.tsv"), "--obo", self.obo])
        self.assertEqual(r.returncode, 2)


class PanelAndRules(unittest.TestCase):
    def test_regions_include_known_allele(self):
        tmp = tempfile.mkdtemp(prefix="t1bed_")
        try:
            bed = os.path.join(tmp, "genes.bed")
            rules = t1common.load_rules(RULES)
            genes = sorted(rules["genes"]["tier_a"] + rules["genes"]["tier_b"])
            with open(bed, "w") as fh:
                for i, g in enumerate(genes):
                    fh.write(f"15\t{1000000 + i * 100000}\t{1050000 + i * 100000}\t{g}\n")
                fh.write("15\t1040000\t1060000\tBUB1B\n")       # 重疊區段要合併
            out = os.path.join(tmp, "regions.bed")
            r = run([os.path.join(SCRIPTS, "panel_regions.py"), "bed", "--bed", bed, "--out", out, "--rules", RULES])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            with open(out) as fh:
                rows = [ln.split("\t") for ln in fh.read().splitlines()]
            self.assertTrue(any(int(s) < 40117088 <= int(e) for c, s, e, _ in rows if c == "15"),
                            "BUB1B 上游已知等位基因要在抽取範圍內")
            starts = [(c, int(s)) for c, s, _, _ in rows]
            self.assertEqual(starts, sorted(starts, key=lambda x: (int(x[0]), x[1])))
            self.assertEqual(len(rows), len(genes) + 1, "重疊區段已合併")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_gene_regex(self):
        r = run([os.path.join(SCRIPTS, "panel_regions.py"), "genes", "--regex", "--rules", RULES])
        names = r.stdout.strip().split("|")
        self.assertIn("CASC5", names)
        self.assertEqual(run([os.path.join(SCRIPTS, "panel_regions.py"), "genes", "--count", "--rules", RULES]).stdout.strip(), "13")

    def test_minrep(self):
        self.assertEqual(t1common.minrep(3000, "AT", "A"), t1common.minrep(3001, "T", "-"))
        self.assertEqual(t1common.minrep(100, "CAG", "CTG"), (101, "A", "T"))

    def test_epcr_ladder(self):
        s = t1common.load_rules(RULES)["submission"]
        lad = s["epcr_primary"] + s["epcr_secondary"]
        self.assertEqual(lad, sorted(set(lad), reverse=True))
        self.assertLessEqual(len(lad), 10)


def bash(code, env):
    e = dict(os.environ)
    e.update(env)
    return subprocess.run(["bash", "-c", f'. "{SCRIPTS}/lib.sh"; {code}'], capture_output=True, text=True, env=e)


class ShellLib(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="t1sh_")
        self.env = {"DATA": os.path.join(self.tmp, "data"), "WORK": os.path.join(self.tmp, "work"),
                    "RESULTS": os.path.join(self.tmp, "results"), "RES": os.path.join(self.tmp, "res")}
        os.makedirs(os.path.join(self.env["DATA"], "sub"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def touch(self, rel):
        p = os.path.join(self.env["DATA"], rel)
        open(p, "w").close()
        return p

    def test_find_vcf(self):
        main = self.touch("sub/WGS_SYNTH.vcf.gz")
        self.touch("sub/WGS_SYNTH.sv.vcf.gz")
        self.touch("sub/cnv_calls.vcf.gz")
        r = bash("find_vcf", self.env)
        self.assertEqual(r.stdout.strip(), main, r.stderr)
        r = bash("sv_vcfs", self.env)
        self.assertEqual(len(r.stdout.split()), 2)
        self.touch("sub/other.vcf.gz")
        self.assertEqual(bash("find_vcf", self.env).stdout.strip(), main, "多個 VCF 時選 WGS_ 開頭的")
        self.touch("sub/WGS_SECOND.vcf.gz")
        r = bash("find_vcf", self.env)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("VCF=", r.stderr)
        e = dict(self.env, VCF=main)
        self.assertEqual(bash("find_vcf", e).stdout.strip(), main)

    def test_fresh(self):
        a, b = os.path.join(self.tmp, "in"), os.path.join(self.tmp, "out")
        for p, t in ((a, "x"), (b, "y")):
            with open(p, "w") as fh:
                fh.write(t)
        os.utime(a, (1000, 1000))
        self.assertEqual(bash(f'fresh "{b}" "{a}"', self.env).returncode, 0)
        os.utime(a, None)
        os.utime(b, (1000, 1000))
        self.assertNotEqual(bash(f'fresh "{b}" "{a}"', self.env).returncode, 0)
        os.utime(b, None)
        self.assertNotEqual(bash(f'fresh "{b}" "{a}"', dict(self.env, FORCE="1")).returncode, 0)

    @unittest.skipUnless(shutil.which("git"), "沒有 git")
    def test_preregistration_gate(self):
        repo = os.path.join(self.tmp, "repo")
        shutil.copytree(T1, os.path.join(repo, "track1"), ignore=shutil.ignore_patterns(".tmp", "__pycache__", "local.env"))
        lib = os.path.join(repo, "track1", "scripts", "lib.sh")
        g = ["git", "-C", repo]
        subprocess.run(g + ["init", "-q"], check=True)
        e = dict(os.environ, **self.env)
        chk = lambda: subprocess.run(["bash", "-c", f'. "{lib}"; check_preregistered'], capture_output=True, text=True, env=e)  # noqa: E731
        r = chk()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("還沒 commit", r.stderr)
        subprocess.run(g + ["add", "track1"], check=True)
        subprocess.run(g + ["-c", "user.name=t", "-c", "user.email=t@example.com", "commit", "-q", "-m", "rules"], check=True)
        self.assertEqual(chk().returncode, 0)
        with open(os.path.join(repo, "track1", "config", "ranking_rules.yml"), "a") as fh:
            fh.write("# changed\n")
        r = chk()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("還沒 commit 的修改", r.stderr)

    @unittest.skipUnless(shutil.which("curl"), "沒有 curl")
    def test_fetch(self):
        """fetch()：伺服器換版不能接成壞檔、連不上不能無聲結束、錯誤頁不能當成檔案大小、.part 要能續傳。"""
        srv_dir = os.path.join(self.tmp, "srv")
        os.makedirs(srv_dir)
        new = b"format-version: 1.2\n" + os.urandom(50000)
        with open(os.path.join(srv_dir, "hp.obo"), "wb") as fh:
            fh.write(new)

        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def send(self, code, body=b"", headers=()):
                self.send_response(code)
                for k, v in headers:
                    self.send_header(k, v)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                if self.command == "GET":
                    self.wfile.write(body)

            def do_HEAD(self):
                self.do_GET()

            def do_GET(self):
                p = os.path.join(srv_dir, self.path.lstrip("/"))
                if not os.path.isfile(p):
                    return self.send(404, b"not found page")
                with open(p, "rb") as fh:
                    data = fh.read()
                m = re.match(r"bytes=(\d+)-$", self.headers.get("Range", ""))
                if m and self.command == "GET":
                    s = int(m.group(1))
                    if s >= len(data):
                        return self.send(416, b"", [("Content-Range", f"bytes */{len(data)}")])
                    return self.send(206, data[s:], [("Content-Range", f"bytes {s}-{len(data) - 1}/{len(data)}")])
                return self.send(200, data)

        try:
            srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        except OSError as e:   # Claude Code 的 sandbox 不讓開本機 port
            self.skipTest(f"沒辦法開本機 port（例如在 Claude Code 的 sandbox 裡）：{e}")
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        self.addCleanup(srv.server_close)
        self.addCleanup(srv.shutdown)
        port = srv.server_address[1]
        sock = socket.socket()
        sock.bind(("127.0.0.1", 0))
        dead = sock.getsockname()[1]          # 沒有人在聽的 port：模擬離線
        sock.close()
        env = {k: v for k, v in os.environ.items() if "proxy" not in k.lower()}
        env.update(self.env, NO_PROXY="127.0.0.1,localhost", no_proxy="127.0.0.1,localhost")

        def fetch(url, out):
            return subprocess.run(["bash", "-c", f'. "{SCRIPTS}/lib.sh"; fetch "$1" "$2"', "_", url, out],
                                  capture_output=True, text=True, env=env)

        def put(path, data):
            with open(path, "wb") as fh:
                fh.write(data)

        def same(path):   # 不用 assertEqual 比 bytes：失敗時 difflib 比對 50 KB 會卡很久
            got = read(path) if os.path.exists(path) else b""
            return got == new, f"得到 {len(got)} bytes，開頭 {got[:8]!r}；應該是 {len(new)} bytes"

        def read(path):
            with open(path, "rb") as fh:
                return fh.read()

        def expect_new(r, path, why=""):
            self.assertEqual(r.returncode, 0, why + "\n" + r.stdout + r.stderr)
            good, detail = same(path)
            self.assertTrue(good, f"{why}：{detail}\n{r.stdout}{r.stderr}")

        url = f"http://127.0.0.1:{port}/hp.obo"
        probe = subprocess.run(["curl", "-s", "-o", os.devnull, "-w", "%{http_code}", "--max-time", "5", url],
                               capture_output=True, text=True, env=env)
        if probe.stdout.strip() != "200":   # sandbox 會擋直接連線到 127.0.0.1
            self.skipTest("連不到本機測試伺服器（例如在 Claude Code 的 sandbox 裡），略過下載測試")
        out = os.path.join(self.env["RES"], "hp.obo")
        os.makedirs(self.env["RES"])
        expect_new(fetch(url, out), out, "第一次下載")
        for old in (b"OLD" + os.urandom(30000), b"OLD" + os.urandom(80000)):   # 伺服器換了新版（舊檔較小、較大）
            put(out, old)
            expect_new(fetch(url, out), out, "換版後必須整個換成新版，不能接成混合檔")
        os.remove(out)
        put(out + ".part", new[:20000])                                        # 上次中斷留下的 .part
        expect_new(fetch(url, out), out, "接續 .part")
        r = fetch(f"http://127.0.0.1:{dead}/hp.obo", out)                      # 檔案完整、伺服器連不上
        expect_new(r, out, "連不上伺服器時要照用已下載的檔案，不能無聲結束")
        self.assertIn("已存在", r.stdout + r.stderr)
        expect_new(fetch(f"http://127.0.0.1:{port}/missing.obo", out), out,   # HEAD 拿到 404 錯誤頁
                   "錯誤頁的 Content-Length 不能拿來核對")
        r = fetch(f"http://127.0.0.1:{dead}/hp.obo", os.path.join(self.env["RES"], "none.obo"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("下載失敗", r.stdout + r.stderr, "沒有檔案又連不上時要印出原因，不能無聲結束")


def mask(text, env=None):
    """用合成字串跑 scripts/mask.sh 的 mva_mask（stdin → stdout）。"""
    e = dict(os.environ)
    if env:
        e.update(env)
    r = subprocess.run(["bash", "-c", f'. "{SCRIPTS}/mask.sh"; mva_mask'],
                       input=text, capture_output=True, text=True, env=e)
    assert r.returncode == 0, r.stderr
    return r.stdout


class Mask(unittest.TestCase):
    """STATUS.md 遮罩（CLAUDE.md 第 5 節）：全部用虛構的合成值，不含真實病人資料。"""

    def test_static_patterns(self):
        s = ("failed at chr15:11112222 and 7:123-456 allele A>G seq ACGTACGT "
             "hgvs c.1234-2A>G p.Ala123Ter g.99998del rs7654321 HP:0009999 "
             "gt 0/1 and 1|0 gene BUB1B alias CASC5 sex MALE count 123456")
        out = mask(s)
        for bad in ["11112222", "A>G", "ACGTACGT", "c.1234", "Ala123", "99998del",
                    "rs7654321", "HP:0009999", "0/1", "1|0", "BUB1B", "CASC5",
                    "MALE", "123456", "chr15:", "7:123"]:
            self.assertNotIn(bad, out, f"沒遮到 {bad!r}：{out}")
        for tag in ["<coord>", "<allele>", "<seq>", "<hgvs>", "<rsid>", "<hpo>",
                    "<gt>", "<gene>", "<sex>", "<num>"]:
            self.assertIn(tag, out, f"少了 {tag}：{out}")

    def test_gene_aliases_and_sex_case(self):
        out = mask("MAD2L1BP FAM178A ZWINT female Female")
        self.assertNotIn("MAD2L1BP", out)
        self.assertNotIn("FAM178A", out)
        self.assertNotIn("ZWINT", out)
        self.assertNotIn("emale", out)   # female / Female 都要遮
        self.assertEqual(out.count("<gene>"), 3, out)
        self.assertEqual(out.count("<sex>"), 2, out)

    def test_runtime_words_and_gated_paths(self):
        out = mask("sample WGS_SYNTH99A at /x/gated-data/sub/WGS_SYNTH99A.vcf.gz done",
                   {"MVA_MASK_WORDS": "WGS_SYNTH99A"})
        self.assertNotIn("WGS_SYNTH99A", out)
        self.assertIn("<masked>", out)
        self.assertIn("gated-data/<masked>", out)

    def test_plain_error_untouched(self):
        s = "bcftools: command not found"
        self.assertEqual(mask(s).strip(), s)
        # 一般句子裡的小數字、時間不用全遮（座標規則以外不動 3 位以下數字）
        self.assertIn("exit code 1", mask("exit code 1"))


class MacCompat(unittest.TestCase):
    """macOS 內建 bash 3.2 與 BSD 指令沒有的寫法，一律不准出現在 .sh 裡。"""
    BAD = [r"\bsed -i\b", r"\bgrep -[a-zA-Z]*P", r"\bstat -c", r"\breadlink -f", r"\bdate -d", r"\bzcat\b",
           r"\bsort -[a-zA-Z]*V", r"\bxargs -r", r"\bdeclare -A", r"\bmapfile\b", r"\breadarray\b", r"\$\{[^}]*,,\}",
           r"\$\{[^}]*\^\^\}", r"&>>", r"\|&", r"-printf\b", r"\btimeout\s", r"\bnproc\b", r"\brealpath\b",
           r"\bhead -n -", r"\btac\b", r"\bmd5sum\b", r"\|\s*grep -[a-zA-Z]*q"]
    # macOS 內建的 bash 3.2 會把緊接在 $var 後面的中文字位元組讀成變數名稱的一部分（set -u 下變成
    # unbound variable）。註解和 heredoc 也檢查：heredoc 會展開變數。
    VAR_THEN_CJK = re.compile(r"\$[A-Za-z_][A-Za-z0-9_]*(?=[^\x00-\x7f])")

    def test_scripts(self):
        for fn in sorted(os.listdir(SCRIPTS)):
            if not fn.endswith(".sh"):
                continue
            path = os.path.join(SCRIPTS, fn)
            r = subprocess.run(["bash", "-n", path], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"{fn}: {r.stderr}")
            raw = open(path, encoding="utf-8").read()
            code = "\n".join(ln.split(" #")[0] for ln in raw.splitlines() if not ln.lstrip().startswith("#"))
            for pat in self.BAD:
                self.assertIsNone(re.search(pat, code), f"{fn} 用了 macOS 沒有的寫法：{pat}")
            m = self.VAR_THEN_CJK.search(raw)
            self.assertIsNone(m, f"{fn}：{m.group(0) if m else ''} 後面直接接中文，macOS 的 bash 3.2 會讀錯變數名稱：改寫成 ${{名稱}}")

    @unittest.skipUnless(shutil.which("shellcheck"), "沒有 shellcheck")
    def test_shellcheck(self):
        files = sorted(os.path.join(SCRIPTS, f) for f in os.listdir(SCRIPTS) if f.endswith(".sh"))
        r = subprocess.run(["shellcheck", "-x", "-S", "warning", "-P", SCRIPTS] + files, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
