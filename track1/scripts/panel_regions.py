#!/usr/bin/env python3
"""MVA 基因組的範圍（基因名單與已知等位基因都來自 config/ranking_rules.yml）。

  panel_regions.py genes --regex | --aliases | --count   給 02_fetch_resources.sh 從 GTF 篩基因
  panel_regions.py bed --bed <基因 BED> --out <輸出>      基因 ±flank + 已知等位基因位置，合併重疊區段
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import default_rules, gene_sets, load_rules, parse_key  # noqa: E402


def chrom_order(c):
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("genes")
    g.add_argument("--rules", default=default_rules())
    m = g.add_mutually_exclusive_group(required=True)
    m.add_argument("--regex", action="store_true", help="現行名稱與舊名，用 | 連接")
    m.add_argument("--aliases", action="store_true", help="舊名=現行名稱，用逗號連接")
    m.add_argument("--count", action="store_true", help="基因數")
    b = sub.add_parser("bed")
    b.add_argument("--rules", default=default_rules())
    b.add_argument("--bed", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--window", type=int, default=50, help="已知等位基因前後各留幾 bp")
    a = ap.parse_args()

    R = load_rules(a.rules)
    tier_a, tier_b, aliases = gene_sets(R)
    genes = sorted(tier_a | tier_b)
    if a.cmd == "genes":
        if a.regex:
            print("|".join(genes + sorted(aliases)))
        elif a.aliases:
            print(",".join(f"{k}={v}" for k, v in sorted(aliases.items())))
        else:
            print(len(genes))
        return

    rows = []
    with open(a.bed, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 4 and not ln.startswith(("#", "track", "browser")):
                rows.append([f[0], int(f[1]), int(f[2]), f[3]])
    found = {r[3] for r in rows}
    if set(genes) - found:
        raise SystemExit(f"{a.bed} 少了基因：{sorted(set(genes) - found)}（重跑 02_fetch_resources.sh panel）")
    for ka in R["known_alleles"]:
        c, p, ref, _ = parse_key(ka["variant"])
        rows.append([c, max(0, p - 1 - a.window), p + len(ref) - 1 + a.window, ka["gene"] + ":known"])
    rows.sort(key=lambda r: (chrom_order(r[0]), r[1], r[2]))
    merged = []
    for r in rows:
        if merged and merged[-1][0] == r[0] and r[1] <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], r[2])
            merged[-1][3] += "," + r[3]
        else:
            merged.append(list(r))
    with open(a.out, "w", encoding="utf-8") as fh:
        for r in merged:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{len(merged)} 個區段（{len(genes)} 個基因 ±flank + {len(R['known_alleles'])} 個已知等位基因）→ {a.out}")


if __name__ == "__main__":
    main()
