#!/usr/bin/env python3
"""Pre-registered sensitivity analysis of route B (config/phase2_protocol.yml, section `sensitivity`).

Re-runs rank_panel.py on a temporary copy of the annotated panel tables, once per scenario, with one rule
changed at a time, and records whether the baseline rank-1 pair stays first. The real route B outputs in
$WORK are never touched, and nothing here changes the submission.

Output (aggregate, no coordinates): JSON with one entry per scenario -> results/track1/report/sensitivity.json
"""
import argparse
import copy
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import POINT_KEYS, default_results, default_rules, default_work, load_rules  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = ["panel.vep.tsv", "panel.spliceai.tsv"]


def read_pairs(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pair_id(row):
    return tuple(sorted(x for x in (row.get("v1", ""), row.get("v2", "")) if x))


def fmt_num(x):
    return f"{x:g}"


def build_scenarios(rules, plan):
    """-> [(group, name, rules_dict)] in a fixed, documented order."""
    out = []
    pts = rules["points"]
    for k in POINT_KEYS:
        for m in plan.get("points_multiply", []):
            r = copy.deepcopy(rules)
            r["points"][k] = pts[k] * m
            out.append(("weights", f"{k} x{fmt_num(m)}", r))
        if plan.get("points_zero"):
            r = copy.deepcopy(rules)
            r["points"][k] = 0
            out.append(("weights", f"{k} removed", r))
    for section, group in (("frequency", "frequency"), ("spliceai", "spliceai"), ("quality", "quality")):
        for key, values in (plan.get(section) or {}).items():
            for v in values:
                r = copy.deepcopy(rules)
                r[section][key] = v
                out.append((group, f"{section}.{key} = {fmt_num(v)}", r))
    for crit in plan.get("tiebreak_first", []):
        r = copy.deepcopy(rules)
        r["pairing"]["tiebreak"] = [crit] + [t for t in rules["pairing"]["tiebreak"] if t != crit]
        out.append(("tiebreak", f"tiebreak: {crit} first", r))
    return out


def run_rank(work_src, rules, tmp, rank_script):
    """Run rank_panel.py on a copy of the inputs with these rules -> (status, pairs or None)."""
    wd = tempfile.mkdtemp(prefix="sens_", dir=tmp)
    for f in INPUTS:
        src = os.path.join(work_src, f)
        if os.path.exists(src):
            shutil.copyfile(src, os.path.join(wd, f))
    rpath = os.path.join(wd, "rules.yml")
    with open(rpath, "w", encoding="utf-8") as fh:
        yaml.safe_dump(rules, fh, allow_unicode=True, sort_keys=False)
    p = subprocess.run([sys.executable, rank_script, "--work", wd, "--rules", rpath, "--top", "1"],
                       capture_output=True, text=True)
    out = os.path.join(wd, "trackB_pairs.tsv")
    if p.returncode != 0 or not os.path.exists(out):
        return "no_pairs", None
    return "ok", read_pairs(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=default_work())
    ap.add_argument("--rules", default=default_rules())
    ap.add_argument("--protocol", default=os.path.join(os.path.dirname(default_rules()), "phase2_protocol.yml"))
    ap.add_argument("--out", default=os.path.join(default_results(), "report", "sensitivity.json"))
    ap.add_argument("--rank-script", default=os.path.join(HERE, "rank_panel.py"))
    a = ap.parse_args()

    rules = load_rules(a.rules)
    with open(a.protocol, encoding="utf-8") as fh:
        proto = yaml.safe_load(fh)
    plan = proto["sensitivity"]["one_at_a_time"]

    base_path = os.path.join(a.work, "trackB_pairs.tsv")
    if not os.path.exists(base_path):
        raise SystemExit(f"Missing {base_path}: run 07_rank.sh first")
    base_rows = read_pairs(base_path)
    if not base_rows:
        raise SystemExit(f"{base_path} is empty")
    base = pair_id(base_rows[0])

    tmp = tempfile.mkdtemp(prefix="t1sens_")
    try:
        status, rows = run_rank(a.work, rules, tmp, a.rank_script)
        reproduced = status == "ok" and pair_id(rows[0]) == base
        results = []
        for group, name, r in build_scenarios(rules, plan):
            status, rows = run_rank(a.work, r, tmp, a.rank_script)
            entry = {"group": group, "name": name}
            if status != "ok":
                entry.update(status="no_pairs", baseline_rank=None, rank1_gene=None)
            else:
                ranks = [i + 1 for i, row in enumerate(rows) if pair_id(row) == base]
                entry.update(status="retained" if ranks and ranks[0] == 1 else "changed",
                             baseline_rank=ranks[0] if ranks else None,
                             rank1_gene=rows[0].get("gene") or None)
            results.append(entry)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    summary = {
        "protocol_version": proto.get("version"),
        "rules_version": rules.get("version"),
        "baseline_reproduced": reproduced,
        "n": len(results),
        "n_retained": sum(1 for x in results if x["status"] == "retained"),
        "n_changed": sum(1 for x in results if x["status"] == "changed"),
        "n_no_pairs": sum(1 for x in results if x["status"] == "no_pairs"),
        "scenarios": results,
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(f"Sensitivity: baseline rank-1 pair kept in {summary['n_retained']}/{summary['n']} scenarios "
          f"(changed {summary['n_changed']}, no pairs {summary['n_no_pairs']}); baseline reproduced: {reproduced}")
    for x in results:
        if x["status"] != "retained":
            print(f"  {x['name']}: {x['status']} (baseline rank {x['baseline_rank']}, rank-1 gene {x['rank1_gene']})")
    print(f"-> {a.out}")
    if not reproduced:
        sys.exit(3)


if __name__ == "__main__":
    main()
