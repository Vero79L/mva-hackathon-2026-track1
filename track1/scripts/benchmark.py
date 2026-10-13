#!/usr/bin/env python3
"""Generalisation benchmark for route B (pre-registered: config/phase2_protocol.yml, section benchmark).

Known pathogenic biallelic genotypes from ClinVar are placed in real public genomes (GIAB HG001-HG007), and
route B, with its frozen code and weights, ranks each genome under eight disease configurations: the MVA rules
and seven PanelApp panels. Public data only; nothing here reads the child's data. 09_benchmark.sh runs the steps:

  format    the sites-only split-vep format, derived from config/splitvep_format.txt
  regions   gene lists + GTF -> one BED per panel, the union BED, panels.json, generated rules per panel
  spikes    ClinVar records (bcftools query output) -> pathogenic variants in each panel's recessive genes
  catalog   split-vep output of the shared annotated catalogue (stdin) -> only rows that could ever qualify
  simulate  backgrounds x panels x genotypes x modes -> rank_panel.py -> benchmark/output/
"""
import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import os
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_rules  # noqa: E402
from rank_panel import COLS, LOF, MISSENSE, num  # noqa: E402  (same parsing as the ranking itself)
from t1common import load_rules, parse_key, strip_chr, variant_key  # noqa: E402

ANN = COLS[13:]                          # SYMBOL ... Existing: what VEP adds, identical for every sample
A_SYM, A_CONS = ANN.index("SYMBOL"), ANN.index("Consequence")
A_MANE, A_CAN = ANN.index("MANE_SELECT"), ANN.index("CANONICAL")
A_AFE, A_AFG, A_MAX = ANN.index("gnomADe_AF"), ANN.index("gnomADg_AF"), ANN.index("MAX_AF")
A_CLN = [ANN.index(c) for c in ("CLNSIG", "CLNREVSTAT", "CLNDN")]
AUTOSOMES = [str(i) for i in range(1, 23)]
MAIN = set(AUTOSOMES) | {"X", "Y"}
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2, "criteria_provided,_single_submitter": 1}
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
HET = dict(GT="0/1", AD="15,15", DP="30", GQ="99")
HOM = dict(GT="1/1", AD="0,30", DP="30", GQ="99")
CLASS_ORDER = {"LoF": 0, "missense": 1, "other": 2}
LABELS = {"mva": "MVA (this study's rules)"}


def opener(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8") if path.endswith(".gz") else open(path, mode, encoding="utf-8")


def bench(protocol):
    with open(protocol, encoding="utf-8") as fh:
        return yaml.safe_load(fh)["benchmark"]


# ---------- format ----------
def sites_format(splitvep_format):
    """Drop the per-sample fields (ID, FILTER, QUAL and the [%...] columns) from the pipeline's split-vep format."""
    with open(splitvep_format, encoding="utf-8") as fh:
        fmt = fh.read().strip()
    toks = [t for t in fmt.replace("\\n", "").split("\\t") if t]
    keep = [t for t in toks if not t.startswith("[") and t not in ("%ID", "%FILTER", "%QUAL")]
    if len(keep) != 4 + len(ANN):
        raise SystemExit(f"{splitvep_format}: expected {4 + len(ANN)} site fields, found {len(keep)}")
    return "\\t".join(keep) + "\\n"


# ---------- regions ----------
def read_gtf_genes(gtf):
    genes = {}
    with opener(gtf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            c = strip_chr(f[0])
            i = f[8].find('gene_name "')
            if c not in MAIN or i < 0:
                continue
            name = f[8][i + 11:f[8].index('"', i + 11)]
            genes.setdefault(name, []).append((c, int(f[3]), int(f[4])))
    return genes


def chrom_order(c):
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def merge_bed(rows):
    rows = sorted(rows, key=lambda r: (chrom_order(r[0]), r[1], r[2]))
    out = []
    for r in rows:
        if out and out[-1][0] == r[0] and r[1] <= out[-1][2]:
            out[-1][2] = max(out[-1][2], r[2])
            out[-1][3] += "," + r[3]
        else:
            out.append(list(r))
    return out


def write_if_changed(path, text):
    """Rewrite only when the content changes, so later steps that compare file times are not redone."""
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            if fh.read() == text:
                return
    with open(path + ".part", "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(path + ".part", path)


def write_bed(rows, path, prefix=""):
    write_if_changed(path, "".join(f"{prefix}{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n" for r in rows))


def cmd_regions(a):
    b = bench(a.protocol)
    base = load_rules(a.rules)
    gtf = read_gtf_genes(a.gtf)
    lists = make_rules.read_genelists(a.genelists)
    names = {}
    if a.meta and os.path.exists(a.meta):
        with open(a.meta, encoding="utf-8") as fh:
            names = {p["key"]: p["name"] for p in json.load(fh)["panels"]}
    os.makedirs(os.path.join(a.out, "rules"), exist_ok=True)
    panels, union = {}, []
    order = ["mva"] + list(b["panels"]["panelapp"])
    for key in order:
        if key == "mva":
            ta, tb = list(base["genes"]["tier_a"]), list(base["genes"]["tier_b"])
            aliases, known = dict(base["genes"].get("aliases") or {}), list(base["known_alleles"])
            biallelic, rules = set(ta), None              # the frozen file itself, not a copy
            src = {"source": "config/ranking_rules.yml", "version": base["version"]}
        else:
            want = b["panels"]["panelapp"][key]
            p = lists.get(key)
            if p is None or str(p["id"]) != str(want["id"]) or p["version"] != str(want["version"]):
                raise SystemExit(f"{a.genelists}: panel {key} must be PanelApp {want['id']} v{want['version']}")
            ta, tb, aliases, known, biallelic = p["tier_a"], p["tier_b"], {}, [], p["biallelic"]
            rules = make_rules.swap(base, ta, tb, note=f"PanelApp {key} (id {p['id']}, version {p['version']})")
            src = {"source": "PanelApp", "id": int(p["id"]), "version": p["version"]}
        rows, missing = [], []
        for g in ta + tb:
            loci = gtf.get(g) or [loc for old, new in aliases.items() if new == g for loc in gtf.get(old, [])]
            if not loci:
                missing.append(g)
            for c, s, e in loci:
                rows.append([c, max(0, s - 1 - a.flank), e + a.flank, g])
        for ka in known:
            c, pos, ref, _ = parse_key(ka["variant"])
            rows.append([c, max(0, pos - 51), pos + len(ref) + 49, ka["gene"] + ":known"])
        spike_genes = sorted(g for g in ta if g in biallelic and any(c in AUTOSOMES for c, _, _ in gtf.get(g, [])))
        bed = merge_bed(rows)
        write_bed(bed, os.path.join(a.out, f"{key}.bed"))
        if rules is None:
            rpath = os.path.abspath(a.rules)
        else:
            rpath = os.path.join(a.out, "rules", f"{key}.yml")
            make_rules.write(rules, rpath + ".new")
            with open(rpath + ".new", encoding="utf-8") as fh:
                write_if_changed(rpath, fh.read())
            os.remove(rpath + ".new")
        union += [r[:3] + [key] for r in bed]
        panels[key] = dict(src, label=LABELS.get(key) or names.get(key, key), tier_a=ta, tier_b=tb, aliases=aliases,
                           known=[ka["variant"] for ka in known], spike_genes=spike_genes, missing_from_gtf=missing,
                           rules=rpath if rules is None else os.path.relpath(rpath, a.out),
                           bp=sum(r[2] - r[1] for r in bed))
        print(f"{key}: tier A {len(ta)}, tier B {len(tb)}, recessive tier-A genes for spikes {len(spike_genes)}, "
              f"{len(bed)} regions, {sum(r[2] - r[1] for r in bed) / 1e6:.1f} Mb"
              + (f"; not in the GTF: {', '.join(missing[:8])}{' …' if len(missing) > 8 else ''}" if missing else ""))
    # backgrounds (GIAB 1-22) and spikes are autosomal, so the shared regions are too
    ubed = merge_bed([r[:4] for r in union if r[0] in AUTOSOMES])
    write_bed(ubed, os.path.join(a.out, "union.bed"))
    write_bed(ubed, os.path.join(a.out, "union.chr.bed"), prefix="chr")
    write_if_changed(os.path.join(a.out, "chrmap.txt"), "".join(f"chr{c}\t{c}\n" for c in AUTOSOMES + ["X", "Y"]))
    write_if_changed(os.path.join(a.out, "panels.json"),
                     json.dumps({"order": order, "panels": panels, "flank": a.flank}, indent=1))
    print(f"shared autosomal regions: {len(ubed)}, {sum(r[2] - r[1] for r in ubed) / 1e6:.1f} Mb -> {a.out}")


# ---------- spikes ----------
def stars(revstat):
    return max([s for t, s in STARS.items() if t in revstat.replace("&", ",")] or [0])


def cmd_spikes(a):
    with open(a.panels, encoding="utf-8") as fh:
        P = json.load(fh)["panels"]
    wanted = {}
    for key, p in P.items():
        for g in p["spike_genes"]:
            wanted.setdefault(g, set()).add(key)
    n_in, kept, seen = 0, [], set()
    with opener(a.raw) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            n_in += 1
            chrom, pos, ref, alt, clnsig, rev, geneinfo, origin = f[:8]
            chrom = strip_chr(chrom)
            sig = clnsig.split("|")[0]
            if chrom not in AUTOSOMES or sig not in PLP or "onflicting" in clnsig or stars(rev) < 1:
                continue
            if max(len(ref), len(alt)) > 50 or set(ref + alt) - set("ACGT") or alt in ("*", "."):
                continue
            o = origin.split(",")[0]
            if o.isdigit() and not int(o) & 1:
                continue                          # ClinVar ORIGIN bit 1 = germline
            genes = sorted({x.split(":")[0] for x in geneinfo.split("|")} & set(wanted))
            k = variant_key(chrom, pos, ref, alt)
            if not genes or k in seen:
                continue
            seen.add(k)
            kept.append((chrom, int(pos), ref, alt, k, ",".join(genes), sig, stars(rev)))
    kept.sort(key=lambda r: (chrom_order(r[0]), r[1], r[2], r[3]))
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("key\tgenes\tclnsig\tstars\n")
        for r in kept:
            fh.write(f"{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\n")
    if a.vcf:
        with open(a.vcf, "w", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            if a.fai:
                with open(a.fai, encoding="utf-8") as fai:
                    for ln in fai:
                        c, length = ln.split("\t")[:2]
                        if c in MAIN:
                            fh.write(f"##contig=<ID={c},length={length}>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for r in kept:
                fh.write(f"{r[0]}\t{r[1]}\t{r[4]}\t{r[2]}\t{r[3]}\t.\t.\t.\n")
    print(f"ClinVar records read {n_in}; eligible pathogenic variants in panel recessive genes {len(kept)} -> {a.out}")


# ---------- catalog ----------
def cmd_catalog(a):
    with open(a.panels, encoding="utf-8") as fh:
        P = json.load(fh)["panels"]
    max_af = load_rules(a.rules)["frequency"]["max_af_recessive"]
    aliases, genes, known = {}, set(), set()
    for p in P.values():
        aliases.update(p["aliases"])
        genes |= set(p["tier_a"]) | set(p["tier_b"])
        known |= {variant_key(*parse_key(k)) for k in p["known"]}
    n_in = n_out = 0
    keys = set()
    src = opener(a.input) if a.input else sys.stdin
    with gzip.open(a.out, "wt", encoding="utf-8") as out:
        for ln in src:
            f = ln.rstrip("\n").split("\t")
            if len(f) != 4 + len(ANN):
                raise SystemExit(f"catalog line has {len(f)} fields, expected {4 + len(ANN)}: use `benchmark.py format`")
            n_in += 1
            ann = f[4:]
            k = variant_key(f[0], f[1], f[2], f[3])
            sym = aliases.get(ann[A_SYM], ann[A_SYM])
            if sym not in genes and k not in known:
                continue
            af = max(num(ann[A_MAX]), num(ann[A_AFE]), num(ann[A_AFG]))
            if af > max_af and "athogenic" not in ann[A_CLN[0]] and k not in known:
                continue                      # common and not ClinVar-pathogenic: can never qualify under the rules
            out.write(ln if ln.endswith("\n") else ln + "\n")
            n_out += 1
            keys.add(k)
    print(f"catalogue rows {n_in} -> kept {n_out} rows for {len(keys)} variants that could qualify -> {a.out}")


# ---------- simulate ----------
def load_catalog(path):
    cat = {}
    with opener(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = variant_key(f[0], f[1], f[2], f[3])
            cat.setdefault(k, []).append(f)
    return cat


def load_calls(path, keys):
    out = {}
    with opener(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            k = variant_key(f[0], f[1], f[2], f[3])
            if k in keys:
                out[k] = f[4:11]                      # ID FILTER QUAL GT AD DP GQ
    return out


def full_row(site, call, pgt=".", pid="."):
    """site = CHROM POS REF ALT + ANN; call = ID FILTER QUAL GT AD DP GQ -> the 30 columns of panel.vep.tsv."""
    return [site[0], site[1], call[0], site[2], site[3], call[1], call[2], call[3], call[4], call[5], call[6],
            pgt, pid] + site[4:]


def preferred(rows, gene, aliases):
    """The transcript row rank_panel.py would keep for (variant, gene): MANE, then canonical, then file order."""
    own = [r for r in rows if aliases.get(r[4 + A_SYM], r[4 + A_SYM]) == gene]
    if not own:
        return None
    return sorted(own, key=lambda r: 0 if len(r[4 + A_MANE]) > 1 else 1 if r[4 + A_CAN] == "YES" else 2)[0]


def vclass(row):
    cons = set(row[4 + A_CONS].split("&"))
    return "LoF" if cons & LOF else "missense" if cons & MISSENSE else "other"


def genotype_class(classes, model):
    c = sorted(classes, key=lambda x: CLASS_ORDER[x])
    return f"hom {c[0]}" if model == "homozygous" else "/".join(c)


def draw(rng, elig, bg_keys, n, p_comp):
    """Seeded draws of spiked genotypes. elig: {gene: [(key, pos, class)]}."""
    free = {g: [v for v in vs if v[0] not in bg_keys] for g, vs in elig.items()}
    comp = sorted(g for g, vs in free.items() if len({v[1] for v in vs}) >= 2)
    hom = sorted(g for g, vs in free.items() if vs)
    out = []
    for _ in range(n):
        want = rng.random() < p_comp
        if want and comp:
            model, pool = "compound_het", comp
        elif hom:
            model, pool = "homozygous", hom
        elif comp:
            model, pool = "compound_het", comp
        else:
            break
        g = rng.choice(pool)
        vs = free[g]
        x = rng.choice(vs)
        if model == "compound_het":
            y = rng.choice([v for v in vs if v[1] != x[1]])
            out.append(dict(model=model, gene=g, keys=[x[0], y[0]], cls=genotype_class([x[2], y[2]], model)))
        else:
            out.append(dict(model=model, gene=g, keys=[x[0]], cls=genotype_class([x[2]], model)))
    return out


def run_task(t):
    """One rank_panel.py run in its own folder. Returns the result row."""
    start = time.time()
    if not t["extra"] and os.path.getsize(t["base"]) == 0:      # nothing in the panel could qualify
        return dict(t["meta"], n_pairs=0, candidates=0, top_points="", rank1_gene="", rank1_points="", seconds=0.0)
    tmp = tempfile.mkdtemp(prefix="bench_", dir=t["tmp"])
    try:
        shutil.copyfile(t["base"], os.path.join(tmp, "panel.vep.tsv"))
        if t["extra"]:
            with open(os.path.join(tmp, "panel.vep.tsv"), "a", encoding="utf-8") as fh:
                fh.write("".join("\t".join(r) + "\n" for r in t["extra"]))
        p = subprocess.run([sys.executable, os.path.join(HERE, "rank_panel.py"), "--work", tmp, "--rules", t["rules"],
                            "--top", "1"], capture_output=True, text=True)
        pairs = []
        path = os.path.join(tmp, "trackB_pairs.tsv")
        if p.returncode == 0 and os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                pairs = list(csv.DictReader(fh, delimiter="\t"))
        elif "沒有合格的配對" not in p.stderr + p.stdout and "沒有任何" not in p.stderr + p.stdout:
            return dict(t["meta"], error=(p.stderr or p.stdout).strip()[-300:], seconds=round(time.time() - start, 2))
        res = dict(t["meta"], n_pairs=len(pairs), seconds=round(time.time() - start, 2),
                   rank1_gene=pairs[0]["gene"] if pairs else "", rank1_points=pairs[0]["pair_points"] if pairs else "")
        cand = [r for r in pairs if r["tier"] == "A" and r["kind"] in ("pair", "hom")]
        res["candidates"] = len(cand)
        res["top_points"] = cand[0]["pair_points"] if cand else ""
        if t["keys"]:
            want = set(t["keys"])
            kind = "hom" if len(want) == 1 else "pair"
            hit = [r for r in pairs if r["kind"] == kind and {r["v1"], r["v2"] or r["v1"]} == want]
            res["rank"] = min(int(r["rank"]) for r in hit) if hit else ""
            res["spike_points"] = hit[0]["pair_points"] if hit else ""
        return res
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def share(xs, pred):
    return round(sum(1 for x in xs if pred(x)) / len(xs), 4) if xs else None


def summary(rows):
    ranks = [r["rank"] for r in rows]
    found = [x for x in ranks if x != ""]
    return dict(n=len(rows), rank1=share(ranks, lambda x: x == 1), top5=share(ranks, lambda x: x != "" and x <= 5),
                top10=share(ranks, lambda x: x != "" and x <= 10), not_found=share(ranks, lambda x: x == ""),
                median_rank=statistics.median(found) if found else None)


def cmd_simulate(a):
    b = bench(a.protocol)
    sp = b["spikes"]
    seed, n_per, p_comp = sp["seed"], int(sp["per_panel_per_background"]), float(sp["genotype_model"]["compound_het"])
    for need, got in (("0/1", HET["GT"]), ("15,15", HET["AD"]), ("1/1", HOM["GT"]), ("0,30", HOM["AD"])):
        if got != need or need not in json.dumps(sp["call_fields"]):
            raise SystemExit(f"protocol call_fields no longer match the code ({need})")
    with open(os.path.join(a.regions, "panels.json"), encoding="utf-8") as fh:
        PJ = json.load(fh)
    P, order = PJ["panels"], PJ["order"]
    max_af = load_rules(a.rules)["frequency"]["max_af_recessive"]
    cat = load_catalog(a.catalog)
    samples = sorted(f[:-7] for f in os.listdir(a.genotypes) if f.endswith(".gt.tsv"))
    if samples != sorted(b["backgrounds"]["samples"]):
        raise SystemExit(f"backgrounds found {samples}, protocol lists {b['backgrounds']['samples']}")
    calls = {s: load_calls(os.path.join(a.genotypes, s + ".gt.tsv"), cat.keys()) for s in samples}
    spikes = {}
    with open(a.spikes, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            spikes[r["key"]] = r
    os.makedirs(a.out, exist_ok=True)
    tmp_root = a.tmp or tempfile.mkdtemp(prefix="bench_")
    os.makedirs(tmp_root, exist_ok=True)
    tasks, panel_info = [], {}
    for key in order:
        p = P[key]
        aliases = p["aliases"]
        genes = set(p["tier_a"]) | set(p["tier_b"])
        known = {variant_key(*parse_key(k)) for k in p["known"]}
        rules = os.path.join(a.regions, p["rules"])

        def relevant(site, genes=genes, aliases=aliases, known=known):
            k = variant_key(site[0], site[1], site[2], site[3])
            return aliases.get(site[4 + A_SYM], site[4 + A_SYM]) in genes or k in known

        # eligible spikes: ClinVar-eligible, annotated on the gene, rare in gnomAD
        elig = {}
        for k, r in spikes.items():
            for g in r["genes"].split(","):
                if g not in p["spike_genes"] or k not in cat:
                    continue
                row = preferred(cat[k], g, aliases)
                if row is None:
                    continue
                af = max(max(num(x[4 + A_MAX]), num(x[4 + A_AFE]), num(x[4 + A_AFG])) for x in cat[k])
                if af > max_af:
                    continue
                c, pos = k.split(":")[:2]
                elig.setdefault(g, []).append((k, (c, int(pos)), vclass(row)))
        for g in elig:
            elig[g].sort()
        panel_info[key] = dict(n_spike_genes_with_variants=len(elig), n_eligible_variants=sum(len(v) for v in elig.values()))
        for s in samples:
            base = os.path.join(tmp_root, f"base_{key}_{s}.tsv")
            with open(base, "w", encoding="utf-8") as fh:
                for k, call in calls[s].items():
                    for site in cat[k]:
                        if relevant(site):
                            fh.write("\t".join(full_row(site, call)) + "\n")
            meta = dict(panel=key, background=s)
            tasks.append(dict(base=base, extra=[], keys=[], rules=rules, tmp=tmp_root, meta=dict(meta, mode="healthy")))
            rng = random.Random(f"{seed}|{key}|{s}")
            for i, gt in enumerate(draw(rng, elig, set(calls[s]), n_per, p_comp)):
                fields = HOM if gt["model"] == "homozygous" else HET
                for mode in ("clinvar", "novel"):
                    extra = []
                    for k in gt["keys"]:
                        c, pos, ref, alt = k.split(":")
                        call = [k, "PASS", "50", fields["GT"], fields["AD"], fields["DP"], fields["GQ"]]
                        for site in cat[k]:
                            if relevant(site):
                                row = full_row(site, call)
                                if mode == "novel":
                                    for j in A_CLN:
                                        row[13 + j] = "."
                                extra.append(row)
                    tasks.append(dict(base=base, extra=extra, keys=gt["keys"], rules=rules, tmp=tmp_root,
                                      meta=dict(meta, mode=mode, i=i, model=gt["model"], gene=gt["gene"],
                                                cls=gt["cls"], spiked=",".join(gt["keys"]))))
    print(f"{len(tasks)} rank_panel.py runs over {len(order)} panels and {len(samples)} genomes, {a.workers} at a time")
    t0 = time.time()
    results = []
    with cf.ProcessPoolExecutor(max_workers=a.workers) as ex:
        for n, res in enumerate(ex.map(run_task, tasks, chunksize=4), 1):
            results.append(res)
            if n % 100 == 0 or n == len(tasks):
                print(f"  {n}/{len(tasks)} ({time.time() - t0:.0f} s)")
    wall = time.time() - t0
    if not a.tmp:
        shutil.rmtree(tmp_root, ignore_errors=True)
    errors = [r for r in results if r.get("error")]
    if errors:
        raise SystemExit(f"{len(errors)} runs failed, first: {errors[0]}")

    cols = ["panel", "background", "mode", "i", "model", "gene", "cls", "spiked", "rank", "spike_points", "n_pairs",
            "candidates", "top_points", "rank1_gene", "rank1_points", "seconds"]
    with open(os.path.join(a.out, "simulations.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in results:
            w.writerow({c: r.get(c, "") for c in cols})

    def sel(**kv):
        return [r for r in results if all(r.get(k) == v for k, v in kv.items())]

    out_panels = []
    for key in order:
        p = P[key]
        healthy = sel(panel=key, mode="healthy")
        entry = dict(key=key, label=p["label"], source=p["source"], id=p.get("id"), version=p.get("version"),
                     n_tier_a=len(p["tier_a"]), n_tier_b=len(p["tier_b"]), n_spike_genes=len(p["spike_genes"]),
                     missing_from_gtf=len(p["missing_from_gtf"]), megabases=round(p["bp"] / 1e6, 2), **panel_info[key],
                     modes={m: summary(sel(panel=key, mode=m)) for m in ("novel", "clinvar")},
                     classes={m: {c: summary([r for r in sel(panel=key, mode=m) if r["cls"] == c])
                                  for c in sorted({r["cls"] for r in sel(panel=key, mode=m)})} for m in ("novel", "clinvar")},
                     healthy=dict(per_background={r["background"]: dict(candidates=r["candidates"], top_points=r["top_points"])
                                                  for r in healthy},
                                  median_candidates=statistics.median([r["candidates"] for r in healthy]) if healthy else None,
                                  max_candidates=max([r["candidates"] for r in healthy]) if healthy else None),
                     median_seconds=statistics.median([r["seconds"] for r in sel(panel=key)]) if sel(panel=key) else None)
        out_panels.append(entry)
        nv = entry["modes"]["novel"]
        print(f"{key:24s} tier A {entry['n_tier_a']:5d}  novel rank1 {nv['rank1']}  top5 {nv['top5']}  "
              f"clinvar rank1 {entry['modes']['clinvar']['rank1']}  healthy candidates median {entry['healthy']['median_candidates']}")
    doc = dict(protocol_version=yaml.safe_load(open(a.protocol, encoding="utf-8")).get("version"),
               protocol_commit=a.protocol_commit, rules_version=load_rules(a.rules)["version"], seed=seed,
               generated=time.strftime("%Y-%m-%d %H:%M"), backgrounds=samples,
               versions=dict(giab="NISTv4.2.1 GRCh38", clinvar=a.clinvar_version, vep=a.vep_version,
                             panelapp={k: P[k].get("version") for k in order if k != "mva"}),
               panels=out_panels,
               overall={m: summary(sel(mode=m)) for m in ("novel", "clinvar")},
               runtime=dict(runs=len(results), wall_seconds=round(wall, 1), workers=a.workers,
                            median_seconds_per_run=statistics.median([r["seconds"] for r in results]),
                            annotation_minutes=a.annotation_minutes if a.annotation_minutes else None))  # 0 = not timed
    with open(os.path.join(a.out, "benchmark.json"), "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
    print(f"-> {os.path.join(a.out, 'benchmark.json')} ({wall / 60:.1f} min)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("format")
    f.add_argument("--splitvep", required=True)
    r = sub.add_parser("regions")
    r.add_argument("--gtf", required=True)
    r.add_argument("--genelists", required=True)
    r.add_argument("--meta")
    r.add_argument("--rules", required=True)
    r.add_argument("--protocol", required=True)
    r.add_argument("--flank", type=int, default=10000)
    r.add_argument("--out", required=True)
    s = sub.add_parser("spikes")
    s.add_argument("--raw", required=True)
    s.add_argument("--panels", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--vcf")
    s.add_argument("--fai")
    c = sub.add_parser("catalog")
    c.add_argument("--panels", required=True)
    c.add_argument("--rules", required=True)
    c.add_argument("--input", help="default: stdin")
    c.add_argument("--out", required=True)
    m = sub.add_parser("simulate")
    m.add_argument("--regions", required=True)
    m.add_argument("--catalog", required=True)
    m.add_argument("--genotypes", required=True)
    m.add_argument("--spikes", required=True)
    m.add_argument("--rules", required=True)
    m.add_argument("--protocol", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    m.add_argument("--tmp")
    m.add_argument("--protocol-commit", default="")
    m.add_argument("--clinvar-version", default="")
    m.add_argument("--vep-version", default="")
    m.add_argument("--annotation-minutes", type=float, default=None)
    a = ap.parse_args()
    if a.cmd == "format":
        print(sites_format(a.splitvep))
    else:
        {"regions": cmd_regions, "spikes": cmd_spikes, "catalog": cmd_catalog, "simulate": cmd_simulate}[a.cmd](a)


if __name__ == "__main__":
    main()
