#!/usr/bin/env python3
"""Fill the Track 1 report template with this run's results.

Runs on the participant's own computer (it reads patient data; it uses no network). It collects every
case-specific value the report needs from git, the run manifest, the pipeline outputs and the two
Phase 2 analyses, writes them as one JSON block into the template, and saves

  results/track1/report/<team>_track1_report.html
  results/track1/report/<team>_track1_report.pdf    (when Chrome, Chromium, Edge or Brave is installed)

Sample names and file paths are never copied into the report.
"""
import argparse
import csv
import datetime as dt
import glob
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
T1 = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(T1, "scripts"))
from t1common import default_results, default_work  # noqa: E402

PANEL_ORDER = ["BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "BUB1", "BUB3", "MAD1L1", "MAD2L1BP",
               "CEP192", "KNL1", "ZWINT"]
ALLELE_POS = {"15:40209701", "15:40220612"}   # GRCh38 positions of c.2210 and c.3006 (NM_001211.6)


# ---------- small readers ----------
def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_text(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return fh.read()
    except OSError:
        return None


def read_json(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def git(t1, *args):
    try:
        p = subprocess.run(["git", "-C", t1] + list(args), capture_output=True, text=True, timeout=30)
        return p.stdout.strip() if p.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def fnum(x):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return None


def to_int(x):
    v = fnum(x)
    return int(v) if v is not None else None


def clean(x):
    return None if x in (None, "", ".", "nan") else x


# ---------- sections ----------
def count_tests(t1):
    n = 0
    tdir = os.path.join(t1, "tests")
    for f in sorted(os.listdir(tdir)) if os.path.isdir(tdir) else []:
        if f.startswith("test_") and f.endswith(".py"):
            n += len(re.findall(r"^\s+def test_", open(os.path.join(tdir, f), encoding="utf-8").read(), re.M))
    return n or None


def prereg(t1, results):
    out = {}
    rc = git(t1, "log", "-1", "--format=%h|%ci", "--", "config/ranking_rules.yml")
    if rc:
        out["rules_commit"], out["rules_commit_time"] = rc.split("|", 1)
    rules_path = os.path.join(t1, "config", "ranking_rules.yml")
    if os.path.exists(rules_path):
        out["rules_sha256"] = hashlib.sha256(open(rules_path, "rb").read()).hexdigest()
        out["rules_version"] = (yaml.safe_load(open(rules_path, encoding="utf-8")) or {}).get("version")
    pc = git(t1, "log", "-1", "--format=%h|%ci", "--", "config/phase2_protocol.yml")
    if pc:
        out["protocol_commit"], out["protocol_time"] = pc.split("|", 1)
    head = git(t1, "log", "-1", "--format=%h|%ci")
    if head:
        out["head_commit"], out["head_time"] = head.split("|", 1)
    # The live copies are in results/; the committed copies in docs/preregistration/ are the fallback.
    def record(name):
        for d in (results, os.path.join(t1, "docs", "preregistration")):
            if os.path.exists(os.path.join(d, name)):
                return os.path.join(d, name)
        return os.path.join(results, name)
    first = read_text(record("first_rank_run.txt"))
    if first:
        m = re.search(r"first ranking run (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})", first)
        if m:
            out["first_rank_run"] = m.group(1)
    out["manifest"] = [{"time": r.get("time"), "step": r.get("step"), "rules_commit": (r.get("rules_commit") or "").split(" ")[0],
                        "note": r.get("note")} for r in read_tsv(record("run_manifest.tsv"))]
    if out.get("first_rank_run"):
        log = git(t1, "log", f"--since={out['first_rank_run']}", "--reverse", "--format=@@%h|%ci|%s", "--name-only", "--", ".")
        changes = []
        for block in (log or "").split("@@")[1:]:
            lines = [x for x in block.strip().split("\n") if x.strip()]
            h, t, subj = lines[0].split("|", 2)
            files = [os.path.relpath(f, "track1") if f.startswith("track1/") else f for f in lines[1:]]
            changes.append({"hash": h, "time": t, "subject": subj, "files": files})
        out["changes_after_first_run"] = changes
    return out


def versions(t1, hpo_obo):
    env = os.path.join(t1, "env")
    v = {"vep": None, "exomiser": None, "exomiser_data": None}
    tools = read_text(os.path.join(env, "tools.txt")) or ""
    for key, pat in (("bcftools", r"^bcftools (\S+)"), ("samtools", r"^samtools (\S+)"), ("python", r"^Python (\S+)"),
                     ("pandas", r"pandas (\S+)"), ("vep", r"ensembl\s+:\s+(\S+)")):
        m = re.search(pat, tools, re.M)
        if m:
            v[key] = m.group(1)
    ex = read_text(os.path.join(env, "exomiser.txt")) or ""
    m = re.search(r"exomiser-cli (\S+), data (\S+)", ex)
    if m:
        v["exomiser"], v["exomiser_data"] = m.group(1).rstrip(","), m.group(2)
    cv = read_text(os.path.join(env, "clinvar_version.txt")) or ""
    m = re.search(r"fileDate=(\S+)", cv)
    if m:
        v["clinvar_date"] = m.group(1)
    jv = read_text(os.path.join(env, "java.txt")) or ""
    m = re.search(r'version "([^"]+)"', jv)
    if m:
        v["java"] = m.group(1)
    mac = read_text(os.path.join(env, "macos.txt")) or ""
    m = re.search(r"ProductVersion:\s*(\S+)", mac)
    if m:
        v["os"] = f"macOS {m.group(1)}"
    sp = read_text(os.path.join(env, "spliceai.yml")) or ""
    m = re.search(r"spliceai[=]{1,2}(\S+)", sp)
    if m:
        v["spliceai"] = m.group(1)
    if hpo_obo and os.path.exists(hpo_obo):
        with open(hpo_obo, encoding="utf-8") as fh:
            for _ in range(30):
                ln = fh.readline()
                if ln.startswith("data-version:"):
                    v["hpo_obo"] = ln.split(":", 1)[1].strip()
                    break
    return {k: x for k, x in v.items() if x}


QC_KEYS = ("records", "snps", "indels", "multiallelic_sites", "ts_tv", "het_hom_ratio", "mean_depth", "pass_percent")


def qc(t1, work):
    """Genome-wide aggregates only. Sex-chromosome metrics stay out of the report, which is published under CC BY."""
    out = {}
    for r in (read_text(os.path.join(t1, "qc", "qc_summary.txt")) or "").splitlines():
        if r.startswith("#") or "\t" not in r:
            continue
        k, x = r.split("\t", 1)
        if k.strip() in QC_KEYS:
            out[k.strip()] = x.strip()
    out["sex_checked"] = bool((read_text(os.path.join(work, "sex.txt")) or "").strip())
    return out


def vcf_info(work):
    t = read_text(os.path.join(work, "vcf_header_check.txt")) or ""
    out = {}
    m = re.search(r"^FORMAT: (.*)$", t, re.M)
    if m:
        fmt = m.group(1).split()
        out["format_fields"] = fmt
        out["has_pgt_pid"] = "PGT" in fmt and "PID" in fmt
    m = re.search(r"^FILTER: (.*)$", t, re.M)
    if m:
        out["filters"] = m.group(1).split()
    return out


def hpo(results):
    """Only the number of terms: the terms themselves are part of the controlled dataset."""
    return {"n": sum(1 for r in read_tsv(os.path.join(results, "hpo_terms.tsv")) if r.get("hpo_id"))}


def exomiser(work, proband):
    presets = []
    for name in ("exome", "introns"):
        path = os.path.join(work, "exomiser", f"{proband}_{name}.genes.tsv")
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        if not rows:
            continue
        head = [h.lstrip("#") for h in rows[0]]
        ix = {h: i for i, h in enumerate(head)}

        def col(r, *names):
            for n in names:
                if n in ix and ix[n] < len(r):
                    return r[ix[n]]
            return None
        genes, seen = [], set()
        for r in rows[1:]:
            g = col(r, "GENE_SYMBOL")
            if not g:
                continue
            genes.append({"rank": to_int(col(r, "RANK")), "gene": g, "moi": col(r, "MOI"),
                          "combined": fnum(col(r, "EXOMISER_GENE_COMBINED_SCORE")),
                          "pheno": fnum(col(r, "EXOMISER_GENE_PHENO_SCORE")),
                          "variant": fnum(col(r, "EXOMISER_GENE_VARIANT_SCORE"))})
        top = []
        for g in genes:                       # one row per gene (best MOI), in Exomiser's order
            if g["gene"] in seen:
                continue
            seen.add(g["gene"])
            top.append(g)
        bub = next((i + 1 for i, g in enumerate(top) if g["gene"] == "BUB1B"), None)
        for i, g in enumerate(top):
            g["rank"] = i + 1
        presets.append({"name": name, "n_genes": len(top), "bub1b_rank": bub, "top": top[:10]})
    return {"presets": presets}


def route_b(work, rules, submitted):
    """Pairs are listed only when every allele is in the submitted CSV (which is public anyway); the rest are counted."""
    pairs = read_tsv(os.path.join(work, "trackB_pairs.tsv"))
    variants = read_tsv(os.path.join(work, "trackB_variants.tsv"))
    excluded = read_tsv(os.path.join(work, "trackB_excluded.tsv"))
    keep = ["rank", "gene", "tier", "kind", "phase", "v1", "v2", "pair_points", "sum_points", "evidence1", "evidence2", "hgvs1", "hgvs2"]
    listed = [p for p in pairs if p.get("v1") in submitted and (not p.get("v2") or p.get("v2") in submitted)]
    cand_genes = sorted({p["gene"] for p in pairs if p.get("kind") in ("pair", "hom")})
    out = {"pairs_top": [{k: p.get(k) for k in keep} for p in listed[:10]], "n_pairs": len(pairs),
           "n_pairs_not_listed": len(pairs) - len(listed), "n_excluded_cis": len(excluded),
           "n_candidate_genes": len(cand_genes), "n_variants": len({v["key"] for v in variants}),
           "n_qualifying": sum(1 for v in variants if v.get("qualifying") == "True")}
    return out, pairs, variants


def alleles(pairs, variants, work):
    if not pairs:
        return []
    top = pairs[0]
    ad = {}
    for r in read_tsv_noheader(os.path.join(work, "panel.vep.tsv")):
        key = f"{r[0]}:{r[1]}:{r[3]}:{r[4]}"
        ad.setdefault(key, r[8] if len(r) > 8 else None)
    out = []
    keys = [k for k in (top.get("v1"), top.get("v2")) if k]
    keys.sort(key=lambda k: (k.split(":")[0], int(k.split(":")[1]) if k.split(":")[1].isdigit() else 0))
    for k in keys:
        if not k:
            continue
        v = next((x for x in variants if x.get("key") == k and x.get("gene") == top.get("gene")), None)
        if not v:
            continue
        out.append({
            "key": k, "gene": v.get("gene"), "hgvsc": clean(v.get("HGVSc")), "hgvsp": clean(v.get("HGVSp")),
            "consequence": clean(v.get("Consequence")), "impact": clean(v.get("IMPACT")), "gt": clean(v.get("GT")),
            "pgt": clean(v.get("PGT")), "pid": v.get("PID") if v.get("PID") not in (None, "") else ".",
            "filter": clean(v.get("FILTER")), "gq": clean(v.get("GQ")), "dp": clean(v.get("DP")), "ad": clean(ad.get(k)),
            "af": fnum(v.get("af")), "spliceai": fnum(v.get("sai")), "clnsig": v.get("CLNSIG") or ".",
            "clnrevstat": clean(v.get("CLNREVSTAT")), "points": fnum(v.get("points")), "evidence": clean(v.get("evidence")),
            "sift": clean(v.get("SIFT")), "polyphen": clean(v.get("PolyPhen")), "transcript": clean(v.get("Feature")),
            "mane": clean(v.get("MANE_SELECT")), "existing": clean(v.get("Existing")),
        })
    return out


def read_tsv_noheader(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]


def panel_genes(variants, work, rules, pairs):
    tier_a = set(rules["genes"]["tier_a"])
    depth = {}
    seen = set()
    aliases = rules["genes"].get("aliases", {})
    for r in read_tsv_noheader(os.path.join(work, "panel.vep.tsv")):
        if len(r) < 14:
            continue
        gene = aliases.get(r[13], r[13])
        key = (gene, f"{r[0]}:{r[1]}:{r[3]}:{r[4]}")
        if gene not in PANEL_ORDER or key in seen:
            continue
        seen.add(key)
        d = fnum(r[9])
        if d is not None:
            depth.setdefault(gene, []).append(d)
    out = []
    cand = {p["gene"] for p in pairs if p.get("kind") in ("pair", "hom")}
    for g in PANEL_ORDER:
        vs = [v for v in variants if v.get("gene") == g]
        ds = depth.get(g, [])
        out.append({"gene": g, "tier": "A" if g in tier_a else "B", "n_variants": len({v["key"] for v in vs}),
                    "candidate": g in cand,
                    "median_dp": statistics.median(ds) if ds else None, "min_dp": min(ds) if ds else None})
    known = next((v for v in variants if v.get("known") == "True"), None)
    return out, ({"present": bool(known and known.get("GT") not in ("0/0", "./.", "0|0")), "gt": known.get("GT") if known else None})


def sv_summary(work):
    rows = read_tsv_noheader(os.path.join(work, "panel_sv.tsv"))
    bed = [r for r in read_tsv_noheader(os.path.join(work, "panel_regions.bed")) if len(r) >= 4]
    items = []
    for r in rows:
        if len(r) < 9:
            continue
        chrom, pos, info = r[1].replace("chr", ""), int(r[2]), r[7]
        m = re.search(r"END=(\d+)", info)
        end = int(m.group(1)) if m else pos + max(len(r[4]), 1)
        t = re.search(r"SVTYPE=([A-Z]+)", info)
        genes = sorted({b[3] for b in bed if b[0].replace("chr", "") == chrom and int(b[1]) < end and int(b[2]) > pos})
        covers = any(f"{chrom}:{p}" in ALLELE_POS and pos <= int(p) <= end for p in ("40209701", "40220612"))
        items.append({"type": t.group(1) if t else "other", "genes": genes, "covers_alleles": covers})
    by_type = {}
    for it in items:
        by_type[it["type"]] = by_type.get(it["type"], 0) + 1
    return {"n": len(items), "by_type": by_type, "genes": sorted({g for it in items for g in it["genes"]}),
            "covers_alleles": any(it["covers_alleles"] for it in items)}


def submission(results, csv_path):
    if not csv_path:
        cands = sorted(glob.glob(os.path.join(results, "*_mva-tiered-v*.csv")), key=os.path.getmtime)
        csv_path = cands[-1] if cands else None
    if not csv_path or not os.path.exists(csv_path):
        return None, {"rows": []}
    rows = []
    with open(csv_path, encoding="utf-8") as fh:
        for i, r in enumerate(csv.DictReader(fh)):
            notes = r.get("notes") or ""
            gene = notes.split(";")[0].split(" (")[0].strip() or None
            row = {k: r.get(k) for k in ("chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type")}
            row.update(row=i + 1, gene=gene, route="A" if "route A" in notes else "B",
                       v1=f"{r['chrom_1']}:{r['pos_1']}:{r['ref_1']}:{r['alt_1']}" if r.get("chrom_1") else None,
                       v2=f"{r['chrom_2']}:{r['pos_2']}:{r['ref_2']}:{r['alt_2']}" if r.get("chrom_2") else None)
            rows.append(row)
    return csv_path, {"rows": rows, "n_secondary": sum(1 for r in rows if (r.get("finding_type") or "").strip() == "secondary")}



# ---------- run time, hardware, cost ----------
STEP_NAMES = {"03_fetch_data": "data download", "04_qc": "QC and sex check", "05_annotate": "normalisation and annotation",
              "06_exomiser": "Exomiser (two presets)", "07_rank": "route B ranking", "08_submission": "CSV and self-check"}
SHORT_NAMES = {"03_fetch_data": "the download", "04_qc": "QC", "05_annotate": "annotation", "06_exomiser": "Exomiser",
               "07_rank": "ranking", "08_submission": "CSV assembly"}
SETUP_STEP = "02_fetch_resources"
LOG_NAME = re.compile(r"(0[2-8]_[a-z_]+?)_(\d{8}-\d{6})\.log$")


def log_failed(path):
    """A run that stopped with die() prints a line starting with the cross mark (lib.sh)."""
    try:
        with open(path, "rb") as fh:
            return any(ln.lstrip().startswith("✗".encode()) for ln in fh)
    except OSError:
        return True


def step_minutes(work):
    """Longest successful run of each pipeline step. Start: the log file's creation time where the file system keeps it
    (macOS), else the time in its name (local time, as lib.sh writes it). End: the log's last write. Runs that stopped
    with an error are skipped. Returns ({step: minutes} for steps 03-08, minutes of the one-time setup step or None)."""
    out = {}
    for path in glob.glob(os.path.join(work, "logs", "*.log")):
        m = LOG_NAME.match(os.path.basename(path))
        if not m or (m.group(1) not in STEP_NAMES and m.group(1) != SETUP_STEP) or log_failed(path):
            continue
        st = os.stat(path)
        start = getattr(st, "st_birthtime", None)
        if not start:
            try:
                start = dt.datetime.strptime(m.group(2), "%Y%m%d-%H%M%S").timestamp()
            except ValueError:
                continue
        mins = max(0.0, (st.st_mtime - start) / 60)
        out[m.group(1)] = max(out.get(m.group(1), 0.0), mins)
    setup = out.pop(SETUP_STEP, None)
    return {k: round(out[k], 1) for k in sorted(out)}, (round(setup, 1) if setup is not None else None)


def sysctl(name):
    try:
        p = subprocess.run(["sysctl", "-n", name], capture_output=True, text=True, timeout=5)
        return p.stdout.strip() if p.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def as_int(x):
    return int(x) if x and x.isdigit() else None


def hardware():
    cpu = sysctl("machdep.cpu.brand_string")
    mem = sysctl("hw.memsize")
    out = {"cpu": cpu, "memory_gb": round(int(mem) / 2 ** 30) if mem and mem.isdigit() else None,
           "cores": as_int(sysctl("hw.physicalcpu")), "perf_cores": as_int(sysctl("hw.perflevel0.physicalcpu")), "os": None}
    mac = sysctl("kern.osproductversion")
    out["os"] = f"macOS {mac}" if mac else (os.uname().sysname if hasattr(os, "uname") else None)
    if not cpu and os.path.exists("/proc/cpuinfo"):
        txt = read_text("/proc/cpuinfo") or ""
        m = re.search(r"model name\s*:\s*(.+)", txt)
        out["cpu"] = m.group(1).strip() if m else None
        out["cores"] = txt.count("processor\t:") or None
    if out.get("perf_cores") == out.get("cores"):
        out["perf_cores"] = None
    return {k: v for k, v in out.items() if v}


def hardware_text(hw):
    cores = None
    if hw.get("cores"):
        cores = f"{hw['cores']} cores" + (f" ({hw['perf_cores']} performance)" if hw.get("perf_cores") else "")
    return ", ".join(x for x in (hw.get("cpu"), cores, f"{hw['memory_gb']} GB memory" if hw.get("memory_gb") else None,
                                 hw.get("os")) if x)


def folder_gb(path):
    """Size of the public resources folder. None when it is missing or nearly empty (not a real resources folder)."""
    if not path or not os.path.isdir(path):
        return None
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            if not os.path.islink(fp):
                try:
                    total += os.path.getsize(fp)
                except OSError:
                    pass
    return round(total / 1e9) if total >= 1e9 else None


def csv_unchanged(csv_path, manifest_rows):
    """True when the CSV was last written by 08_submission.sh (its time matches the run manifest), else False."""
    if not csv_path or not os.path.exists(csv_path):
        return None
    times = []
    for r in manifest_rows:
        if (r.get("step") or "").startswith("08") and (r.get("note") or "") in ("", os.path.basename(csv_path)):
            try:
                times.append(dt.datetime.strptime(r["time"], "%Y-%m-%d %H:%M:%S").timestamp())
            except (KeyError, ValueError, TypeError):
                pass
    if not times:
        return None
    return os.path.getmtime(csv_path) <= max(times) + 120


# ---------- the organisers' methods form ----------
FORM_TITLE = "Methods description — MVA Hackathon 2026, Track 1 (Variant Prediction)"
FORM_INTRO = ("Please submit your answers with your Track 1 submission. All questions (with the exception of Question 1, in red) "
              "are voluntary; however, submissions without a methods description may be judged less favorably on the "
              "Innovation and Scalability criteria, since the panel can't assess what it can't see.\n\n"
              "Fill out this form once per model/approach submitted.")
FORM_Q = [
    (7, "Team name"),
    (8, "Model number (fill out this form once per model/approach submitted; up to 6)"),
    (9, "Please describe your model/approach in detail."),
    (10, "(required) If commercially-available Generative AI was used, record the provider, the plan or tier, and the "
         "relevant setting in your methods description. A line is enough, for example, \"Anthropic API, Claude Sonnet, "
         "commercial terms, no training on customer content.\""),
    (11, "Please state whether the submission file is the automated output of your computational approach (preferred), "
         "or if it has undergone downstream manual review and curation."),
    (12, "Please describe any downstream manual review in detail."),
    (13, "Please state whether your approach only used publicly available data (preferred) or if proprietary data was "
         "also used. If proprietary data were used, where possible, please also include a submission based only on "
         "publicly available data."),
    (14, "Please describe any public data used in detail (e.g., reference panels, population databases, literature)."),
    (15, "Please describe any proprietary data used in detail."),
    (16, "Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single "
         "candidate variants?"),
    (17, "How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)"),
    (18, "Please provide an estimate of run time and cost for your approach (if possible)."),
    (19, "Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well "
         "as methodological detail."),
]


def pct(x):
    return f"{round(100 * x)}%" if isinstance(x, (int, float)) else "n/a"


def cls_word(a):
    cons = (a or {}).get("consequence") or ""
    if any(c in cons for c in ("stop_gained", "frameshift", "splice_donor_variant", "splice_acceptor_variant", "start_lost")):
        return "truncating"
    if any(c in cons for c in ("missense", "inframe")):
        return "missense"
    return "non-coding" if cons else "unclassified"


def pooled_rank1(panels):
    """Share ranked first across panels, weighted by the number of runs (as in the report's text)."""
    n = sum(p["modes"]["novel"]["n"] for p in panels)
    return sum(p["modes"]["novel"]["rank1"] * p["modes"]["novel"]["n"] for p in panels) / n if n else None


def bench_sentence(b):
    if not b:
        return None
    pan = {p["key"]: p for p in b["panels"]}
    tested = [p for p in b["panels"] if p["modes"]["novel"]["n"]]
    others = [p for p in tested if p["key"] != "mva"]
    mva = pan.get("mva", {}).get("modes", {}).get("novel", {})
    sizes = [p["n_tier_a"] for p in others]
    small = [p for p in tested if p["n_tier_a"] <= 30]
    large = [p for p in tested if p["n_tier_a"] > 100]
    by_size = (f" By panel size: first in {pct(pooled_rank1(small))} for panels of up to 30 green genes and in "
               f"{pct(pooled_rank1(large))} for panels above 100." if small and large else "")
    return (f"A pre-registered benchmark placed ClinVar-pathogenic biallelic genotypes into {len(b['backgrounds'])} public "
            f"genomes (GIAB) and re-ran route B, unchanged except for the gene list, under the MVA rules and "
            f"{len(others)} PanelApp disease panels ({min(sizes):,} to {max(sizes):,} green genes). " if sizes else
            f"A pre-registered benchmark placed ClinVar-pathogenic biallelic genotypes into {len(b['backgrounds'])} public "
            f"genomes (GIAB) and re-ran route B under the MVA rules. ") + (
            f"With ClinVar evidence removed, as for newly found variants, the placed genotype ranked first "
            f"in {pct(b['overall']['novel']['rank1'])} of {b['overall']['novel']['n']} runs (MVA rules: {pct(mva.get('rank1'))}) "
            f"and in the top five in {pct(b['overall']['novel']['top5'])}; with ClinVar evidence, first in "
            f"{pct(b['overall']['clinvar']['rank1'])}.{by_size}")


def named(label, val):
    return f"{label} {val}" if val else label


def nice_date(x):
    return f"{x[:4]}-{x[4:6]}-{x[6:8]}" if x and re.fullmatch(r"\d{8}", x) else x


def minutes_text(x):
    return "under 1 min" if x < 0.5 else f"{round(x)} min"


def duration_text(mins):
    return f"{round(mins)} minutes" if mins < 90 else f"{mins / 60:.1f} hours"


def article(word):
    return "an" if word[:1] in "aeiou" else "a"


def plural(n, one, many):
    return f"{'no' if not n else n} {one if n == 1 else many}"


def rules_unchanged(pre):
    """True when the ranking rules were last committed before the first ranking run; None when that is unknown."""
    t, first = pre.get("rules_commit_time"), pre.get("first_rank_run")
    if not t or not first:
        return None
    try:
        tc = dt.datetime.strptime(t, "%Y-%m-%d %H:%M:%S %z").astimezone().replace(tzinfo=None)
        return tc <= dt.datetime.strptime(first, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        return None


def written_for_this_pair(al):
    """The report's evidence narrative, and the abstract's classification sentence, were written for these two alleles."""
    p = " ".join((x.get("hgvsp") or "") for x in al)
    return len(al) == 2 and bool(re.search(r"Leu737(Ter|\*)", p)) and "Asn1002Lys" in p


def form_answers(d, a):
    m, v = d["meta"], d["versions"]
    sub, rb, ra = d["submission"], d["routeB"], d["routeA"]
    al = d.get("alleles") or []
    sens, an, b = d.get("sensitivity"), d.get("aneuploidy"), d.get("benchmark")
    rt, hw = d.get("runtime") or {}, d.get("hardware") or {}
    train = m.get("ai_training") or "<on or off>"
    n_genes = len(d.get("panel_genes") or []) or 13
    n_tests = m.get("n_tests")
    presets = ra.get("presets") or []
    ranks = [p.get("bub1b_rank") for p in presets if p.get("bub1b_rank")]
    top_gene = ((rb.get("pairs_top") or [{}])[0].get("gene") if rb.get("pairs_top") else None) or "BUB1B"
    route_a = ("Route A (Exomiser) also ranked BUB1B first in both presets." if ranks and all(r == 1 for r in ranks) and len(ranks) == 2
               else f"Route A ranked BUB1B {', '.join(str(r) for r in ranks)} in its presets." if ranks else "")
    classes = sorted(cls_word(x) for x in al) if len(al) == 2 else []
    rows = sub.get("rows") or []
    r1 = rows[0] if rows else {}
    row1 = ("The submitted row 1 is a compound-heterozygous pair." if r1.get("v1") and r1.get("v2") and r1["v1"] != r1["v2"]
            else "The submitted row 1 is a homozygous variant." if r1.get("v1") and r1.get("v2") == r1.get("v1")
            else "The submitted row 1 is a single variant." if r1.get("v1") else "")
    steps = rt.get("steps") or {}
    total = sum(steps.values()) if steps else None
    hw_txt = hardware_text(hw) or "one laptop"
    rg, setup_min = d.get("resources_gb"), rt.get("setup_minutes")
    unchanged = d.get("csv_unchanged")
    n_sec = sub.get("n_secondary", 0)
    frozen = rules_unchanged(d.get("prereg") or {})
    b_sent = bench_sentence(b)
    sens_txt = (f"A pre-registered sensitivity analysis changed one rule at a time in {sens.get('n')} scenarios; the "
                f"rank-1 pair stayed first in {sens.get('n_retained')}." if sens and sens.get("n") else "")
    an_txt = ("A pre-registered scan of the child's genome for mosaic whole-chromosome aneuploidy was "
              + ("negative" if not any(x.get("positive") for x in an.get("prespecified", [])) and not an.get("exploratory_flags") else "reported")
              + ", with its detection limit stated." if an else "")
    vep, clinvar = named("Ensembl VEP", v.get("vep")), "ClinVar" + (f" release {nice_date(v['clinvar_date'])}" if v.get("clinvar_date") else "")
    exo = named("Exomiser", v.get("exomiser")) + (f" with data release {v['exomiser_data']}" if v.get("exomiser_data") else "")
    sai = named("SpliceAI", v.get("spliceai"))
    ans = {
        7: f"{m.get('team') or '<team>'} (Hugging Face username; individual participant).",
        8: f"{m.get('model') or '1'}",
        9: (f"Two independent routes share normalisation ({named('bcftools', v.get('bcftools'))}) and annotation ({vep} with "
            f"MANE Select transcripts and gnomAD v4.1 frequencies, {clinvar}, {sai}) and then rank separately. Route A: {exo} "
            f"ranks every gene in the genome from the child's {plural(d['hpo'].get('n'), 'HPO term', 'HPO terms')}, with the "
            "exome and exome-plus-introns presets, and makes no assumption about the diagnosis. Route B scores every variant, "
            f"coding and non-coding, PASS and non-PASS, in {n_genes} MVA and spindle-checkpoint genes with a points system "
            "adapted from the Bayesian ACMG/AMP framework, adds a known upstream hypomorph, pairs alleles within each gene, "
            "scores a pair by its weaker allele and excludes pairs proven to be in cis by read-backed phasing. Its weights, "
            "thresholds and tie-breaks were committed to git before the first ranking run"
            + (" and have not changed since." if frozen else "; later changes are listed from git in the report." if frozen is False
               else "; the report lists the commit history.")
            + " Code merges both routes into the 10-row CSV with a fixed EPCR ladder and checks it with the official scorer. "
            "Phase 2 adds per-allele ACMG/AMP evidence, a pre-registered sensitivity analysis, a mosaic-aneuploidy scan of the "
            "child's genome and a pre-registered benchmark on public genomes in which only the gene list changes between "
            "diseases. Full detail is in the report sections Design, Pipeline, Robustness, Aneuploidy and Scalability."),
        10: (f"Anthropic Claude (claude.ai and Claude Code), Max subscription, consumer terms, model training setting {train}. "
             "Used to write code, tests and documentation and to draft the report. No AI tool opened the data folders "
             "(Claude Code ran with deny rules and an operating-system sandbox; claude.ai was connected to code folders "
             "only), and every step that reads data ran as a plain script started by the participant. Candidate-level "
             "results were pasted into conversation after the analysis had finished; this is disclosed in the report "
             "(AI use and data governance). Claude Agent SDK and the Anthropic API were not used."),
        11: ("Automated output of make_submission.py, which merges route B and route A by fixed rules and assigns fixed "
             "EPCR values. " + ("The file has not been changed since the code wrote it: its modification time matches "
                                 "the submission step in the run manifest." if unchanged else
                                 "The file was changed after the code wrote it; see the next answer." if unchanged is False else "")),
        12: ("No row was added, removed, reordered or edited by hand. Rows 9 and 10 are reserved for secondary findings, "
             "which must be reviewed by a person and listed with a reason in a curated file before the code adds them; "
             f"{plural(n_sec, 'such row was', 'such rows were')} submitted. The run book also asks for a look at the "
             "rank-1 pair's call metrics and the structural-variant overlap before upload; that review changes nothing "
             "in the file."),
        13: "Only publicly available data and open-source software, together with the organisers' dataset for this child "
            "(the VCF and the HPO terms). No proprietary data.",
        14: (f"GRCh38 reference (Ensembl); {vep} cache with MANE Select transcripts, gnomAD v4.1 allele frequencies, SIFT "
             f"and PolyPhen-2; {clinvar}; {sai} with its published models; {exo} (HPO annotations, disease models, gnomAD, "
             "ClinVar, protein interactions); Human Phenotype Ontology"
             + (f" (hp.obo {v['hpo_obo']})" if v.get("hpo_obo") else "") + "; UniProt; about 70 published papers, "
             "guidelines and database records cited in the report (ACMG/AMP and ClinGen calibrations, MVA literature, SIOP-Europe surveillance "
             "guidance). Benchmark only: GIAB NISTv4.2.1 genomes HG001 to HG007 and Genomics England PanelApp panels at "
             "pinned versions."),
        15: "None.",
        16: ("Yes. Route B outputs biallelic candidates directly: compound-heterozygous pairs within a gene (pairs proven "
             "to be in cis are excluded; unphased pairs are kept and labelled), homozygous variants, and, ranked after "
             "them, genes with a single qualifying allele. Route A adds gene-level candidates under recessive models. "
             + row1),
        17: ("Secondary findings are kept apart from the primary answer: rows 9 and 10 of the CSV are reserved for them, "
             "they need a person's review and a written reason in a curated file, and they are flagged "
             f"finding_type = secondary. {plural(n_sec, 'secondary row was', 'secondary rows were')} submitted. Incidental "
             "results are not reported in the public write-up."),
        18: ((f"Measured on {hw_txt}. " if steps or hw else "Estimate for one laptop. ") + "One-time setup: download of "
             + (f"about {rg} GB of " if rg else "the ") + "public resources listed in answer 14"
             + (f", which took {duration_text(setup_min)} on this connection" if setup_min else "") + ". Per genome: "
             + ("; ".join(f"{STEP_NAMES[k]} {minutes_text(x)}" for k, x in steps.items())
                + f"; about {duration_text(total)} in all, most of it {' and '.join(SHORT_NAMES[k] for k in sorted(steps, key=steps.get, reverse=True)[:2])}"
                if steps else "about 1.5 to 3 hours on a laptop, most of it Exomiser (estimate; step logs not available)")
             + ". The Phase 2 analyses and the report add a few minutes. "
             + (f"In the benchmark, annotating all sites took {max(1, round(b['runtime']['annotation_minutes']))} minutes once, and "
                if b and (b.get("runtime") or {}).get("annotation_minutes") else "In the benchmark, " if b else "")
             + (f"each route B ranking took a median of {b['runtime']['median_seconds_per_run']} s ({b['runtime']['runs']} runs). "
                if b else "")
             + "No cloud computing or paid service; every tool and dataset is free and open, so the marginal cost per "
               "genome is laptop time and electricity. Genomes are processed independently, so a cohort scales linearly "
               "and parallelises across machines. No AI model is part of the pipeline."),
    }
    if len(al) == 2:
        first, second = (classes[1], classes[0]) if classes[0] == "missense" else (classes[0], classes[1])
        result = (f"Result. Route B ranked first a {top_gene} genotype made of {article(first)} {first} allele and "
                  f"{article(second)} {second} allele"
                  + (", the configuration of most published MVA1 families" if top_gene == "BUB1B" and classes == ["missense", "truncating"] else "")
                  + f". {route_a} "
                  + ("Under ClinGen-refined ACMG/AMP criteria the truncating allele classifies as pathogenic; the missense "
                     "allele remains of uncertain significance until parental testing confirms trans phase and a protein "
                     "assay is done." if written_for_this_pair(al) else "Per-allele ACMG/AMP evidence is set out in the report."))
    else:
        result = f"Result. Route B's top candidate is set out in the report. {route_a}"
    scale = ("Scalability. The ranking engine and the disease knowledge are separate files: another recessive disease needs "
             "a new gene list, which one command builds from a PanelApp panel, not new code. " + (b_sent or "")
             + (f" One genome takes about {duration_text(total)} on a laptop, without cloud or paid services, and genomes "
                "run independently, so a cohort scales linearly." if total else " Genomes run independently on a laptop, "
                "so a cohort scales linearly."))
    abstract = "\n\n".join(" ".join(x.split()) for x in [
        "Approach. We ranked candidate variants for a child with mosaic variegated aneuploidy (MVA) along two routes "
        "that share normalisation and annotation but not assumptions. Route A runs Exomiser genome-wide from the "
        f"child's phenotype terms. Route B scores every variant in {n_genes} MVA and spindle-checkpoint genes with a points "
        "system adapted from the Bayesian ACMG/AMP framework, keeps non-coding and non-PASS calls, adds a known "
        "upstream hypomorph, pairs alleles within each gene, scores a pair by its weaker allele and excludes pairs "
        "proven to be in cis. Every weight, threshold and tie-break was committed to git before the first ranking run, "
        "and every later commit is listed from git. Code merges the routes into the submission and checks it with the "
        "official scorer.",
        result,
        f"Robustness. {sens_txt} {an_txt}" if (sens_txt or an_txt) else "",
        scale,
        "Strengths. Pre-registration with an audit trail; two independent routes; recessive-aware pairing with cis "
        "exclusion; non-coding and non-PASS variants in scope; per-allele evidence on calibrated scales; "
        + ("generalisation tested on public genomes with only the gene list changed; " if b else "")
        + "laptop-scale run time" + (f"; {n_tests} automated tests." if n_tests else "."),
        "Limitations. Phase is inferred, not observed, so trans confirmation needs the parents. Route B depends on a "
        "correct gene list, and with large panels background candidates compete, which is why it is paired with the "
        "gene-agnostic route A. " + ("The benchmark uses public calls with synthetic call metrics and no SpliceAI. " if b else "")
        + "Mobile-element insertions, distant regulatory variants and low-level mosaic variants are largely invisible in "
        "a short-read germline VCF.",
    ] if x.strip())
    ans[19] = abstract
    return [{"row": r, "q": qtext, "a": " ".join(ans[r].split()) if r != 19 else ans[r].strip(), "required": r == 10}
            for r, qtext in FORM_Q], len(abstract.split())


def xml_text(t):
    return (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;"))


def write_form_xlsx(form, path):
    """The organisers' form layout (same rows) with column B filled in. Standard library only."""
    import zipfile
    widths = {"A": 60, "B": 140}
    rows = {1: [("A", FORM_TITLE, 1)], 2: [("A", FORM_INTRO, 0)], 6: [("B", "Answers", 1)]}
    for f in form:
        rows[f["row"]] = [("A", f["q"], 2 if f["required"] else 1), ("B", f["a"].replace("\n\n", "\n"), 0)]

    def height(cells):
        """Row height in points for wrapped text, so every app shows the whole answer (Excel caps rows at 409 pt)."""
        lines = max(sum(-(-max(len(par), 1) // widths[col]) for par in txt.split("\n")) for col, txt, _ in cells)
        return min(409, 15 * lines + 3)
    sheet = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
             '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
             f'<cols><col min="1" max="1" width="{widths["A"]}" customWidth="1"/>'
             f'<col min="2" max="2" width="{widths["B"]}" customWidth="1"/></cols>',
             "<sheetData>"]
    for r in sorted(rows):
        sheet.append(f'<row r="{r}" ht="{height(rows[r])}" customHeight="1">' + "".join(
            f'<c r="{col}{r}" t="inlineStr" s="{st}"><is><t xml:space="preserve">{xml_text(txt)}</t></is></c>'
            for col, txt, st in rows[r]) + "</row>")
    sheet.append("</sheetData></worksheet>")
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              '<fonts count="3"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font>'
              '<font><b/><sz val="11"/><color rgb="FFC00000"/><name val="Calibri"/></font></fonts>'
              '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
              '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
              '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
              '<cellXfs count="3">'
              '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf>'
              '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf>'
              '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf>'
              '</cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>')
    files = {
        "[Content_Types].xml": '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>',
        "_rels/.rels": '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        "xl/workbook.xml": '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Track 1 methods" sheetId="1" r:id="rId1"/></sheets></workbook>',
        "xl/_rels/workbook.xml.rels": '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>',
        "xl/worksheets/sheet1.xml": "".join(sheet),
        "xl/styles.xml": styles,
    }
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, body in files.items():
            z.writestr(name, body)


def write_form_md(form, n_words, path):
    lines = [f"# {FORM_TITLE}", "", "Answers in the organisers' form, generated by track1/report/fill_report.py. "
             "The answers describe the method; they contain no genotype-level data.", ""]
    for f in form:
        lines += [f"## {f['row']}. {f['q']}", "", f["a"], ""]
    lines.append(f"_Abstract length: {n_words} words._")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

# ---------- output ----------
def write_html(template, data, out_html):
    html = open(template, encoding="utf-8").read()
    payload = json.dumps(data, ensure_ascii=False, indent=1).replace("</", "<\\/")
    new, n = re.subn(r'(<script id="report-data" type="application/json">)(.*?)(</script>)',
                     lambda m: m.group(1) + payload + m.group(3), html, count=1, flags=re.S)
    if n != 1:
        raise SystemExit("The template has no report-data block")
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(new)


def find_browser():
    env = os.environ.get("CHROME")
    cands = [env] if env else []
    cands += ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "/Applications/Chromium.app/Contents/MacOS/Chromium",
              "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge", "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser"]
    cands += [shutil.which(n) for n in ("google-chrome", "chromium", "chromium-browser", "chrome")]
    return next((c for c in cands if c and os.path.exists(c)), None)


def write_pdf(out_html, out_pdf):
    browser = find_browser()
    if not browser:
        print("No Chrome/Chromium/Edge/Brave found: open the HTML and use Print > Save as PDF.")
        return False
    prof = tempfile.mkdtemp(prefix="t1pdf_")
    try:
        cmd = [browser, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check", f"--user-data-dir={prof}",
               "--no-pdf-header-footer", "--print-to-pdf-no-header", "--virtual-time-budget=15000", "--run-all-compositor-stages-before-draw",
               f"--print-to-pdf={out_pdf}", Path(out_html).resolve().as_uri() + "#print"]
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            cmd.insert(1, "--no-sandbox")      # only needed when run as root (containers)
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        ok = p.returncode == 0 and os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 10000
        if not ok:
            print(f"PDF export failed ({os.path.basename(browser)}): {p.stderr.strip()[-300:]}")
        return ok
    finally:
        shutil.rmtree(prof, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=default_work())
    ap.add_argument("--results", default=default_results())
    ap.add_argument("--t1", default=T1)
    ap.add_argument("--template", default=os.path.join(HERE, "report_template.html"))
    ap.add_argument("--out-dir", default=None, help="default: <results>/report")
    ap.add_argument("--csv", default=None, help="submission CSV (default: newest *_mva-tiered-v*.csv in results)")
    ap.add_argument("--proband", default=os.environ.get("PROBAND", "PROBAND01"))
    ap.add_argument("--hpo-obo", default=os.environ.get("HPO_OBO"))
    ap.add_argument("--team", default=None, help="default: the CSV file name prefix (HF username)")
    ap.add_argument("--author", default=None)
    ap.add_argument("--github", default=None)
    ap.add_argument("--ai-training", default=None, help='claude.ai "Help improve Claude" setting: on or off')
    ap.add_argument("--lb-submission", default=None)
    ap.add_argument("--lb-rank-points", default=None)
    ap.add_argument("--lb-fmax", default=None)
    ap.add_argument("--lb-match", default=None, help="full / partial / none")
    ap.add_argument("--model", default="1", help="model number in the organisers' form (one form per model)")
    ap.add_argument("--res", default=os.environ.get("RES"), help="public resources folder, for the size estimate")
    ap.add_argument("--no-form-in-repo", action="store_true", help="do not copy the filled methods form into track1/docs")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args()

    with open(os.path.join(a.t1, "config", "ranking_rules.yml"), encoding="utf-8") as fh:
        rules = yaml.safe_load(fh)
    csv_path, sub = submission(a.results, a.csv)
    team = a.team or (os.path.basename(csv_path).split("_mva-tiered")[0] if csv_path else None)
    submitted = {k for r in sub["rows"] for k in (r.get("v1"), r.get("v2")) if k}
    submitted |= {k.replace("chr", "", 1) for k in submitted}
    rb, pairs, variants = route_b(a.work, rules, submitted)
    genes, known = panel_genes(variants, a.work, rules, pairs)
    out_dir = a.out_dir or os.path.join(a.results, "report")
    an = read_json(os.path.join(out_dir, "aneuploidy.json"))
    if an:                                   # autosomes only in the report; sex chromosomes stay out of the public text
        an = {k: v for k, v in an.items() if k != "sex"}
        an["chroms"] = [c for c in an.get("chroms", []) if c.get("chrom") not in ("X", "Y")]
    pre = prereg(a.t1, a.results)
    steps, setup_min = step_minutes(a.work)
    data = {
        "meta": {"team": team, "author": a.author, "github": a.github, "ai_training": a.ai_training, "model": a.model,
                 "generated": dt.datetime.now().strftime("%Y-%m-%d %H:%M"), "csv_name": os.path.basename(csv_path) if csv_path else None,
                 "n_tests": count_tests(a.t1)},
        "prereg": pre,
        "versions": versions(a.t1, a.hpo_obo),
        "vcf": vcf_info(a.work),
        "qc": qc(a.t1, a.work),
        "hpo": hpo(a.results),
        "routeA": exomiser(a.work, a.proband),
        "routeB": rb,
        "alleles": alleles(pairs, variants, a.work),
        "panel_genes": genes,
        "known_allele": known,
        "sv": sv_summary(a.work),
        "submission": sub,
        "sensitivity": read_json(os.path.join(out_dir, "sensitivity.json")),
        "aneuploidy": an,
        "benchmark": read_json(os.path.join(a.t1, "benchmark", "output", "benchmark.json")),
        "runtime": {"steps": steps, "setup_minutes": setup_min},
        "hardware": hardware(),
        "resources_gb": folder_gb(a.res),
        "csv_unchanged": csv_unchanged(csv_path, pre.get("manifest", [])),
        "leaderboard": {"submission": a.lb_submission, "rank_points": a.lb_rank_points, "f_max": a.lb_fmax, "match": a.lb_match},
        "rules_text": read_text(os.path.join(a.t1, "config", "ranking_rules.yml")),
        "protocol_text": read_text(os.path.join(a.t1, "config", "phase2_protocol.yml")),
    }
    data["form"], data["form_words"] = form_answers(data, a)
    stem = f"{team or 'team'}_track1_report"
    out_html = os.path.join(out_dir, stem + ".html")
    form_stem = f"{team or 'team'}_track1_methods_form"
    form_dirs = [out_dir] + ([] if a.no_form_in_repo else [os.path.join(a.t1, "docs")])
    for d in form_dirs:
        os.makedirs(d, exist_ok=True)
        write_form_md(data["form"], data["form_words"], os.path.join(d, form_stem + ".md"))
        write_form_xlsx(data["form"], os.path.join(d, form_stem + ".xlsx"))
    print(f"Methods form: {', '.join(os.path.join(d, form_stem + '.xlsx') for d in form_dirs)} ({data['form_words']} words in the abstract)")
    if data["form_words"] > 500:
        print("WARNING: the abstract is over the organisers' 500-word limit; shorten it in fill_report.py (form_answers).")
    if len(data["alleles"]) == 2 and not written_for_this_pair(data["alleles"]):
        print("WARNING: the report's evidence narrative was written for p.(Leu737Ter) and p.(Asn1002Lys); "
              "the top pair differs, so revise section 'Variant evidence' before submitting.")
    write_html(a.template, data, out_html)
    with open(os.path.join(out_dir, stem + ".data.json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
    print(f"Report: {out_html}")
    missing = [k for k in ("alleles", "routeA", "sensitivity", "aneuploidy", "benchmark") if not data[k] or data[k] == {"presets": []}]
    if missing:
        print("Not yet available (shown as placeholders): " + ", ".join(missing))
    if not a.no_pdf:
        pdf = os.path.join(out_dir, stem + ".pdf")
        if write_pdf(out_html, pdf):
            print(f"PDF:    {pdf}")


if __name__ == "__main__":
    main()
