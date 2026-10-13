#!/usr/bin/env python3
"""Pre-registered scan for mosaic whole-chromosome aneuploidy in the child's own WGS
(config/phase2_protocol.yml, section `aneuploidy`).

Two statistics per chromosome, each compared with the other autosomes by a robust z-score:
  depth     mean allele depth (ref + alt) at PASS SNVs, divided by the autosomal mean (sites outside
            min_depth .. max_depth_factor x autosomal median are dropped first)
  balance   mean z^2 at heterozygous SNVs, z = (alt - n*p0) / sqrt(n*p0*(1-p0))
A mosaic gain or loss raises both; copy-neutral LOH raises balance only.

Memory stays small (compact arrays and histograms), so a 30x genome runs on a laptop in a few minutes.
Output (aggregate per chromosome, no positions): results/track1/report/aneuploidy.json
"""
import argparse
import json
import math
import os
import statistics
import subprocess
import sys
from array import array

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import default_results, default_rules, default_work  # noqa: E402

AUTOSOMES = [str(i) for i in range(1, 23)]
CHROMS = AUTOSOMES + ["X", "Y"]
PAR = [(10001, 2781479), (155701383, 156030895)]   # GRCh38 pseudo-autosomal regions on X
MAXV = 65535
FBINS = 2000                                       # allele-fraction histogram resolution (0.0005)


def in_par(pos):
    return any(a <= pos <= b for a, b in PAR)


def read_sites(vcf, bcftools):
    """Stream PASS biallelic SNVs -> {chrom: (depths, het depths, het alt counts)} as uint16 arrays."""
    view = subprocess.Popen([bcftools, "view", "-f", "PASS", "-v", "snps", "-m2", "-M2", "-Ou", vcf],
                            stdout=subprocess.PIPE)
    query = subprocess.Popen([bcftools, "query", "-f", "%CHROM\t%POS\t%ID\t[%GT]\t[%AD]\n"],
                             stdin=view.stdout, stdout=subprocess.PIPE, text=True)
    view.stdout.close()
    data = {c: (array("H"), array("H"), array("H")) for c in CHROMS}
    for ln in query.stdout:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        chrom = f[0][3:] if f[0].startswith("chr") else f[0]
        if chrom not in data or "," in f[2]:
            continue                       # records split from a multi-allelic site keep the original ALT list in ID
        if chrom == "X" and in_par(int(f[1])):
            continue
        gt = f[3].replace("|", "/")
        ad = f[4].split(",")
        if len(ad) != 2 or gt not in ("0/1", "1/0", "1/1", "1"):
            continue                       # "1": haploid call (some callers on male X/Y); used for depth only
        try:
            r, a = int(ad[0]), int(ad[1])
        except ValueError:
            continue
        n = min(r + a, MAXV)
        dep, hn, ha = data[chrom]
        dep.append(n)
        if gt not in ("1/1", "1"):
            hn.append(n)
            ha.append(min(a, MAXV))
    query.wait()
    view.wait()
    if query.returncode != 0 or view.returncode != 0:
        raise SystemExit(f"bcftools failed on {vcf}")
    return data


def hist_median(hist):
    total = sum(hist)
    if not total:
        return None
    half, run = (total + 1) / 2, 0
    for v, c in enumerate(hist):
        run += c
        if run >= half:
            return v
    return len(hist) - 1


def robust(values):
    med = statistics.median(values)
    mad = statistics.median([abs(v - med) for v in values])
    return med, 1.4826 * mad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", default=os.path.join(default_work(), "norm.vcf.gz"))
    ap.add_argument("--sex", default="UNKNOWN", help="MALE / FEMALE / UNKNOWN (from work/track1/sex.txt)")
    ap.add_argument("--protocol", default=os.path.join(os.path.dirname(default_rules()), "phase2_protocol.yml"))
    ap.add_argument("--out", default=os.path.join(default_results(), "report", "aneuploidy.json"))
    ap.add_argument("--bcftools", default="bcftools")
    a = ap.parse_args()

    with open(a.protocol, encoding="utf-8") as fh:
        full = yaml.safe_load(fh)
    proto = full["aneuploidy"]
    min_dp = proto["sites"]["min_depth"]
    max_factor = proto["sites"]["max_depth_factor"]
    z_pre = proto["prespecified"]["z_threshold"]
    z_exp = proto["exploratory"]["z_threshold"]

    data = read_sites(a.vcf, a.bcftools)
    hist = [0] * (MAXV + 1)
    for c in AUTOSOMES:
        for d in data[c][0]:
            hist[d] += 1
    if sum(hist) < 1000:
        raise SystemExit("Too few autosomal PASS SNVs to run the scan")
    med_depth = hist_median(hist)
    hi = max_factor * med_depth

    def keep(n):
        return min_dp <= n <= hi

    fh_ = [0] * (FBINS + 1)
    for c in AUTOSOMES:
        for n, al in zip(data[c][1], data[c][2]):
            if keep(n):
                fh_[int(al / n * FBINS)] += 1
    p0 = (hist_median(fh_) + 0.5) / FBINS
    q0 = 1 - p0

    auto_sum = auto_n = 0
    for c in AUTOSOMES:
        for d in data[c][0]:
            if keep(d):
                auto_sum += d
                auto_n += 1
    auto_mean = auto_sum / auto_n

    rows = []
    for c in CHROMS:
        dep, hn, ha = data[c]
        n_sites = dsum = 0
        for d in dep:
            if keep(d):
                dsum += d
                n_sites += 1
        z2 = n_het = nsum = 0
        for n, al in zip(hn, ha):
            if keep(n):
                z2 += (al - n * p0) ** 2 / (n * p0 * q0)
                n_het += 1
                nsum += n
        rows.append({
            "chrom": c,
            "n_sites": n_sites,
            "n_het": n_het,
            "depth_ratio": round(dsum / n_sites / auto_mean, 4) if n_sites >= 200 else None,
            "balance": round(z2 / n_het, 4) if n_het >= 200 else None,
            "mean_n": round(nsum / n_het, 2) if n_het else None,
        })

    auto = [r for r in rows if r["chrom"] in AUTOSOMES and r["depth_ratio"] is not None and r["balance"] is not None]
    if len(auto) < 15:
        raise SystemExit("Too few autosomes with enough heterozygous sites to estimate the null spread")
    med_d, sd_d = robust([r["depth_ratio"] for r in auto])
    med_b, sd_b = robust([r["balance"] for r in auto])
    sd_d = max(sd_d, 0.002)     # floors keep z finite when chromosomes agree almost exactly
    sd_b = max(sd_b, 0.002)
    tested = set(AUTOSOMES) | ({"X"} if a.sex == "FEMALE" else set())   # male X/Y: descriptive only
    for r in rows:
        ok = r["chrom"] in tested
        r["z_depth"] = round((r["depth_ratio"] - med_d) / sd_d, 2) if ok and r["depth_ratio"] is not None else None
        r["z_balance"] = round((r["balance"] - med_b) / sd_b, 2) if ok and r["balance"] is not None else None
        flag = None
        if r["chrom"] in AUTOSOMES and r["z_depth"] is not None and r["z_balance"] is not None:
            if r["z_balance"] >= z_exp and r["z_depth"] >= z_exp:
                flag = "gain"
            elif r["z_balance"] >= z_exp and r["z_depth"] <= -z_exp:
                flag = "loss"
            elif r["z_balance"] >= z_exp and abs(r["z_depth"]) < 2:
                flag = "cn_loh"
        r["flag"] = flag

    pre = []
    for c in proto["prespecified"]["chromosomes"]:
        r = next(x for x in rows if x["chrom"] == str(c))
        ok = r["z_depth"] is not None and r["z_balance"] is not None
        pre.append({"chrom": str(c), "z_depth": r["z_depth"], "z_balance": r["z_balance"],
                    "positive": bool(ok and r["z_depth"] >= z_pre and r["z_balance"] >= z_pre)})

    mean_n = statistics.median([r["mean_n"] for r in auto])
    delta = math.sqrt(3 * sd_b / (4 * mean_n))          # allele-fraction shift that moves `balance` by 3 SD
    f_depth = min(1.0, 2 * 3 * sd_d)                      # depth ratio is 1 +/- f/2 for a one-copy change
    power = {
        "sd_depth_ratio": round(sd_d, 5),
        "sd_balance": round(sd_b, 5),
        "mean_het_depth": mean_n,
        "min_cell_fraction_gain": round(min(1.0, max(f_depth, 4 * delta / max(1e-9, 1 - 2 * delta))), 3),
        "min_cell_fraction_loss": round(min(1.0, max(f_depth, 4 * delta / (1 + 2 * delta))), 3),
        "min_cell_fraction_cnloh": round(min(1.0, 2 * delta), 3),
    }
    out = {
        "protocol_version": full.get("version"),
        "sex": a.sex,
        "p0": round(p0, 4),
        "median_depth": med_depth,
        "mean_depth": round(auto_mean, 2),
        "n_sites_autosomal": sum(r["n_sites"] for r in auto),
        "n_het_autosomal": sum(r["n_het"] for r in auto),
        "chroms": rows,
        "prespecified": pre,
        "exploratory_flags": [{"chrom": r["chrom"], "flag": r["flag"]} for r in rows if r["flag"]],
        "power": power,
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"Aneuploidy scan: {out['n_het_autosomal']} autosomal heterozygous SNVs, median depth {med_depth}, p0 {p0:.3f}")
    print("Pre-specified (gain of 7 or 8): " + ", ".join(
        f"chr{x['chrom']} z_depth {x['z_depth']} z_balance {x['z_balance']} -> {'POSITIVE' if x['positive'] else 'negative'}"
        for x in pre))
    print("Exploratory flags: " + (", ".join(f"chr{x['chrom']} {x['flag']}" for x in out["exploratory_flags"]) or "none"))
    print(f"Detectable cell fraction: gain {power['min_cell_fraction_gain']:.0%}, loss {power['min_cell_fraction_loss']:.0%}, "
          f"CN-LOH {power['min_cell_fraction_cnloh']:.0%}")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
