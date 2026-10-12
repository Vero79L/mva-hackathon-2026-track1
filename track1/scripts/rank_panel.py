#!/usr/bin/env python3
"""路線 B：MVA 基因組的事先登錄排序（Track 1 · Phase 1 · Spec v1.2）

輸入（都在 $WORK，含基因型，不進 repo）
  panel.vep.tsv        05_annotate.sh 的 split-vep 輸出（欄位順序 = config/splitvep_format.txt）
  panel.spliceai.tsv   選配；沒有時 SpliceAI 分數以 0 計
輸出
  trackB_variants.tsv  每個（變異 × 基因）一列：點數、證據、是否合格
  trackB_pairs.tsv     同一基因內的配對，依規則排好（rank 欄）
  trackB_excluded.tsv  因為 PGT/PID 證實在同一條染色體（cis）而排除的配對

所有分數、門檻、決勝順序都來自 config/ranking_rules.yml；這支程式裡沒有可調的數字。
"""
import argparse
import csv
import itertools
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import default_rules, default_work, gene_sets, load_rules, parse_key, variant_key  # noqa: E402

COLS = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER", "QUAL", "GT", "AD", "DP", "GQ", "PGT", "PID",
        "SYMBOL", "Feature", "Consequence", "IMPACT", "MANE_SELECT", "CANONICAL", "HGVSc", "HGVSp",
        "SIFT", "PolyPhen", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "CLNSIG", "CLNREVSTAT", "CLNDN", "Existing"]
LOF = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost",
       "transcript_ablation"}
MISSENSE = {"missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "stop_lost"}
NONCODING = {"upstream_gene_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
             "splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant",
             "splice_polypyrimidine_tract_variant", "regulatory_region_variant", "TF_binding_site_variant"}
HET = {"0/1", "1/0", "0|1", "1|0"}
HOM = {"1/1", "1|1"}
IMPACT = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1}


def num(x):
    """VEP 的頻率欄可能是空白、'.'、或用 & 連接的多個值；取最大值，沒有就是 0。"""
    vals = []
    for t in str(x).replace(",", "&").split("&"):
        try:
            vals.append(float(t))
        except ValueError:
            pass
    return max(vals) if vals else 0.0


def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def read_vep(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"找不到或是空檔：{path}（先跑 05_annotate.sh）")
    with open(path, encoding="utf-8") as fh:
        n = fh.readline().rstrip("\n").count("\t") + 1
    if n != len(COLS):
        raise SystemExit(f"{path} 有 {n} 欄，預期 {len(COLS)} 欄：config/splitvep_format.txt 和 rank_panel.py 的 COLS 要一致")
    return pd.read_csv(path, sep="\t", names=COLS, dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)


def read_spliceai(path, aliases):
    """→ {(變異, 基因): 四個 delta 分數的最大值}。只採用和變異同一個 ALT 的註解，基因名稱換成現行名稱。"""
    out = {}
    if not path or not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5 or "|" not in f[4]:
                continue
            key = variant_key(f[0], f[1], f[2], f[3])
            for entry in f[4].split(","):
                p = entry.split("|")
                if len(p) < 6 or p[0].upper() != f[3].upper():
                    continue
                ds = [float(x) for x in p[2:6] if x not in ("", ".")]
                if ds:
                    k = (key, aliases.get(p[1], p[1]))
                    out[k] = max(out.get(k, 0.0), max(ds))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=default_work())
    ap.add_argument("--rules", default=default_rules())
    ap.add_argument("--top", type=int, default=15, help="畫面顯示前幾名配對")
    a = ap.parse_args()

    R = load_rules(a.rules)
    P, F, S, Q = R["points"], R["frequency"], R["spliceai"], R["quality"]
    TIER_A, TIER_B, ALIASES = gene_sets(R)
    PANEL = TIER_A | TIER_B
    KNOWN = {}
    for ka in R["known_alleles"]:
        c, p, r, al = parse_key(ka["variant"])
        KNOWN[variant_key(c, p, r, al)] = ka

    # 1. 讀 VEP 表；基因名稱換成現行名稱
    v = read_vep(os.path.join(a.work, "panel.vep.tsv"))
    v["key"] = [variant_key(c, p, r, x) for c, p, r, x in zip(v.CHROM, v.POS, v.REF, v.ALT)]
    v["gene"] = v.SYMBOL.map(lambda s: ALIASES.get(s, s))
    own = v[v.gene.isin(PANEL)].assign(own=True)

    # 已知等位基因：可能在基因間（VEP 不會標成 MVA 基因），用 rules 指定的基因
    extra = []
    for k, ka in KNOWN.items():
        if ((own.key == k) & (own.gene == ka["gene"])).any():
            continue
        m = v[v.key == k]
        if len(m):
            extra.append(m.iloc[[0]].assign(gene=ka["gene"], own=False, IMPACT="MODIFIER", HGVSc="", HGVSp=""))
    allv = pd.concat([own] + extra, ignore_index=True)
    if allv.empty:
        raise SystemExit("panel.vep.tsv 裡沒有任何 MVA 基因的變異：檢查 05_annotate.sh 的 panel 抽取與 VEP 註解")

    # 2. 每個（變異 × 基因）只留一個轉錄本：MANE > canonical > 其他
    allv["tx"] = 2
    allv.loc[allv.CANONICAL == "YES", "tx"] = 1
    allv.loc[allv.MANE_SELECT.str.len() > 1, "tx"] = 0
    allv = allv.sort_values(["key", "gene", "tx"], kind="stable").drop_duplicates(["key", "gene"]).reset_index(drop=True)

    # 3. SpliceAI、頻率、ClinVar、品質
    sai = read_spliceai(os.path.join(a.work, "panel.spliceai.tsv"), ALIASES)
    if sai is None:
        print("（找不到 panel.spliceai.tsv：SpliceAI 分數以 0 計，第二階段補跑）")
        sai = {}
    allv["sai"] = [sai.get((k, g), 0.0) for k, g in zip(allv.key, allv.gene)]
    allv["af"] = [max(num(m), num(e), num(g)) for m, e, g in zip(allv.MAX_AF, allv.gnomADe_AF, allv.gnomADg_AF)]
    allv["plp"] = allv.CLNSIG.str.contains("athogenic") & ~allv.CLNSIG.str.contains("onflicting")
    allv["stars"] = allv.CLNREVSTAT.map(lambda x: max([s for t, s in STARS.items() if t in x.replace("&", ",")] or [0]))
    allv["known"] = [k in KNOWN and KNOWN[k]["gene"] == g for k, g in zip(allv.key, allv.gene)]
    allv["lowqual"] = [f != "PASS" or (to_int(gq) is not None and to_int(gq) < Q["min_gq"])
                       or (to_int(dp) is not None and to_int(dp) < Q["min_dp"])
                       for f, gq, dp in zip(allv.FILTER, allv.GQ, allv.DP)]
    allv["impact"] = allv.IMPACT.map(lambda x: IMPACT.get(x, 0))

    # 4. 逐一變異打分（只有本基因轉錄本上的後果才算本基因的證據）
    def score(r):
        pts, ev = 0, []
        cons = set(r.Consequence.split("&")) if r.own else set()
        tier_a = r.gene in TIER_A
        rare = r.af <= F["max_af_recessive"]
        if tier_a and cons & LOF:
            pts += P["PVS1"]; ev.append("PVS1")
        cv, cv_tag = 0, ""
        if r.own and r.plp and r.stars >= 2:
            cv, cv_tag = P["clinvar_plp_2star"], "ClinVar_PLP_2star+"
        elif r.own and r.plp and r.stars == 1:
            cv, cv_tag = P["clinvar_plp_1star"], "ClinVar_PLP_1star"
        if r.known and P["known_allele"] > cv:
            cv, cv_tag = P["known_allele"], "Known_allele"
        if cv:
            pts += cv; ev.append(cv_tag)
        if r.sai >= S["strong"]:
            pts += P["spliceai_strong"]; ev.append("SpliceAI_strong")
        elif r.sai >= S["supporting"]:
            pts += P["spliceai_supporting"]; ev.append("SpliceAI_supporting")
        if "missense_variant" in cons and "deleterious" in r.SIFT and "damaging" in r.PolyPhen:
            pts += P["pp3_sift_polyphen"]; ev.append("PP3")
        if tier_a and rare and cons & MISSENSE:
            pts += P["tierA_missense_rare"]; ev.append("TierA_missense_rare")
        if tier_a and rare and cons & NONCODING and not cons & (LOF | MISSENSE):
            pts += P["tierA_noncoding_rare"]; ev.append("TierA_noncoding_rare")
        functional = bool(ev)
        if r.af < F["pm2_max_af"]:
            pts += P["pm2_supporting"]; ev.append("PM2_supporting")
        if tier_a:
            pts += P["pp4_tierA"]; ev.append("PP4")
        return pts, ";".join(ev), functional

    res = [score(r) for r in allv.itertuples(index=False)]
    allv["points"] = [x[0] for x in res]
    allv["evidence"] = [x[1] for x in res]
    allv["functional"] = [x[2] for x in res]
    rare_ok = (allv.af <= F["max_af_recessive"]) | (allv.plp & bool(F["clinvar_plp_bypass"])) | allv.known
    allv["qualifying"] = rare_ok & allv.functional & allv.GT.isin(HET | HOM)
    allv["tier"] = allv.gene.map(lambda g: "A" if g in TIER_A else "B")
    keep = ["key", "gene", "tier", "own", "qualifying", "points", "evidence", "GT", "PGT", "PID", "FILTER", "GQ",
            "DP", "lowqual", "Consequence", "IMPACT", "HGVSc", "HGVSp", "af", "sai", "CLNSIG", "CLNREVSTAT",
            "known", "ID", "Feature", "MANE_SELECT", "SIFT", "PolyPhen", "Existing"]
    allv[keep].sort_values(["tier", "gene", "points"], ascending=[True, True, False], kind="stable") \
        .to_csv(os.path.join(a.work, "trackB_variants.tsv"), sep="\t", index=False)

    # 5. 同一基因內配對
    def allele(r):
        hgvs = next((h for h in (r.HGVSp, r.HGVSc) if h not in ("", ".")), r.key)
        return dict(key=r.key, points=r.points, impact=r.impact, lowqual=r.lowqual, af=r.af,
                    evidence=r.evidence, hgvs=hgvs)

    def make(gene, tier, kind, phase, alle, pair_points):
        w = min(alle, key=lambda x: (x["points"], x["impact"]))
        return dict(gene=gene, tier=tier, kind=kind, phase=phase,
                    v1=alle[0]["key"], v2=alle[1]["key"] if len(alle) > 1 and kind == "pair" else "",
                    pair_points=round(pair_points, 3), weaker_impact=w["impact"],
                    phase_rank={"trans": 1, "hom": 1, "unphased": 0, "cis": -1}.get(phase, 0),
                    sum_points=sum(x["points"] for x in alle),
                    quality_n=sum(0 if x["lowqual"] else 1 for x in alle),
                    max_af=max(x["af"] for x in alle),
                    evidence1=alle[0]["evidence"], evidence2=alle[1]["evidence"] if kind == "pair" else "",
                    hgvs1=alle[0]["hgvs"], hgvs2=alle[1]["hgvs"] if kind == "pair" else "")

    pairs, excluded = [], []
    q = allv[allv.qualifying]
    for gene, g in q.groupby("gene", sort=True):
        tier = "A" if gene in TIER_A else "B"
        homs = [r for r in g.itertuples(index=False) if r.GT in HOM]
        hets = [r for r in g.itertuples(index=False) if r.GT in HET]
        for h in homs:
            al = allele(h)
            pairs.append(make(gene, tier, "hom", "hom", [al, al], al["points"]))
        for x, y in itertools.combinations(hets, 2):
            if x.ID == y.ID and "," in x.ID:
                phase = "trans"                      # 原本是同一個多等位基因位點的 1/2 → 必定 trans
            elif x.PID not in ("", ".") and x.PID == y.PID:
                phase = "cis" if x.PGT == y.PGT else "trans"
            else:
                phase = "unphased"
            ax, ay = allele(x), allele(y)
            row = make(gene, tier, "pair", phase, [ax, ay], min(ax["points"], ay["points"]))
            if phase == "cis" and R["pairing"]["exclude_cis"]:
                excluded.append(row)
                continue
            pairs.append(row)
        if tier == "A" and len(hets) == 1 and not homs:
            al = allele(hets[0])
            pairs.append(make(gene, tier, "single_het_second_allele_missing", "single", [al],
                              al["points"] * R["pairing"]["single_het_factor"]))

    if not pairs:
        raise SystemExit("沒有合格的配對：檢查非 PASS 變異、SV/CNV（05 的 panel_sv.tsv），第二階段用 FASTQ 找第二個等位基因")

    tb = {"weaker_impact": ("weaker_impact", False), "phase_trans": ("phase_rank", False),
          "sum_points": ("sum_points", False), "quality": ("quality_n", False), "rarer": ("max_af", True)}
    out = pd.DataFrame(pairs)
    out["tier_rank"] = (out.tier != "A").astype(int)
    out["kind_rank"] = out.kind.str.startswith("single").astype(int)
    by = ["tier_rank", "kind_rank", "pair_points"] + [tb[t][0] for t in R["pairing"]["tiebreak"]] + ["v1", "v2"]
    asc = [True, True, False] + [tb[t][1] for t in R["pairing"]["tiebreak"]] + [True, True]
    out = out.sort_values(by, ascending=asc, kind="stable").drop(columns=["tier_rank", "kind_rank"])
    out.insert(0, "rank", range(1, len(out) + 1))
    out.to_csv(os.path.join(a.work, "trackB_pairs.tsv"), sep="\t", index=False)
    pd.DataFrame(excluded).to_csv(os.path.join(a.work, "trackB_excluded.tsv"), sep="\t", index=False)

    # 6. 摘要（只印在畫面與 $WORK 的 log，不進 repo）
    print(f"規則 version {R['version']}｜MVA 基因組變異 {allv.key.nunique()} 個｜合格等位基因 {int(allv.qualifying.sum())} 個"
          f"｜配對 {len(out)} 組｜cis 排除 {len(excluded)} 組")
    per = allv.groupby("gene").agg(variants=("key", "size"), qualifying=("qualifying", "sum")).reset_index()
    print("各基因：" + "、".join(f"{r.gene} {int(r.qualifying)}/{r.variants}" for r in per.itertuples()))
    show = ["rank", "gene", "tier", "kind", "phase", "pair_points", "v1", "v2", "hgvs1", "hgvs2"]
    with pd.option_context("display.width", 200, "display.max_colwidth", 40):
        print(out[show].head(a.top).to_string(index=False))
    print(f"\n輸出：{a.work}/trackB_pairs.tsv、trackB_variants.tsv、trackB_excluded.tsv")


if __name__ == "__main__":
    main()
