#!/usr/bin/env python3
"""整理 Exomiser（路線 A）的基因層級結果：每次執行的前 N 名基因、MVA 基因各排第幾。

  exomiser_summary.py --dir <work>/exomiser --out <results>/exomiser_summary.md [--top 10]

只輸出基因名稱、名次、遺傳模式和分數（屬於可保留的「基因排名」），不輸出變異座標。
"""
import argparse
import csv
import datetime
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import default_rules, gene_sets, load_rules  # noqa: E402


def find_col(header, *names):
    norm = {h.lstrip("#").strip().upper(): h for h in header}
    for n in names:
        if n.upper() in norm:
            return norm[n.upper()]
    return None


def read_genes(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return [], {}
    h = list(rows[0].keys())
    cols = {"rank": find_col(h, "RANK"), "gene": find_col(h, "GENE_SYMBOL", "GENE"),
            "moi": find_col(h, "MOI", "MODE_OF_INHERITANCE"),
            "combined": find_col(h, "EXOMISER_GENE_COMBINED_SCORE", "COMBINED_SCORE"),
            "pheno": find_col(h, "EXOMISER_GENE_PHENO_SCORE", "PHENO_SCORE"),
            "variant": find_col(h, "EXOMISER_GENE_VARIANT_SCORE", "VARIANT_SCORE")}
    if not cols["rank"] or not cols["gene"]:
        raise SystemExit(f"{path} 的欄位和預期不同：{h[:12]}")
    return rows, cols


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rules", default=default_rules())
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args()
    tier_a, tier_b, _ = gene_sets(load_rules(a.rules))

    files = sorted(glob.glob(os.path.join(a.dir, "*.genes.tsv")))
    if not files:
        raise SystemExit(f"{a.dir} 裡沒有 *.genes.tsv：先跑 06_exomiser.sh")
    md = ["# Exomiser（路線 A）基因層級摘要", "",
          f"產生時間：{datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}。只含基因名次與分數，不含變異座標。", ""]
    tops = {}
    for f in files:
        run = os.path.basename(f).replace(".genes.tsv", "")
        rows, c = read_genes(f)
        md += [f"## {run}", "", "| 名次 | 基因 | 遺傳模式 | 綜合 | 表現型 | 變異 |", "|---|---|---|---|---|---|"]
        for r in rows[:a.top]:
            md.append("| {} | {} | {} | {} | {} | {} |".format(
                r[c["rank"]], r[c["gene"]], r.get(c["moi"], "") if c["moi"] else "",
                *(("%.3f" % float(r[c[k]])) if c[k] and r.get(c[k]) not in (None, "") else "" for k in ("combined", "pheno", "variant"))))
        best = {}
        for r in rows:
            g = r[c["gene"]]
            if g in tier_a | tier_b and g not in best:
                best[g] = (int(float(r[c["rank"]])), r.get(c["moi"], "") if c["moi"] else "")
        md += ["", "MVA 基因的最佳名次：" + ("、".join(
            f"{g}（{'A' if g in tier_a else 'B'}）#{best[g][0]} {best[g][1]}" for g in sorted(best, key=lambda x: best[x][0]))
            or "前面的結果裡都沒有 MVA 基因"), ""]
        tops[run] = rows[0][c["gene"]] if rows else ""
    md += ["## 兩條路線是否一致", "",
           "各次執行的第 1 名基因：" + "、".join(f"{k} → {v}" for k, v in tops.items()) +
           "。和路線 B（trackB_pairs.tsv）第 1 名的基因相同時，在報告標記「兩條路線一致」。", ""]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
