#!/usr/bin/env python3
"""組提交 CSV、格式檢查、本機自我評分（Track 1 · Spec v1.2 第 11 節）

  make_submission.py draft [--name <HF 帳號>]
      依 ranking_rules.yml 的 merge 順序與 EPCR 階梯，從路線 B（trackB_pairs.tsv）、路線 A（Exomiser）
      和人工審核的 secondary findings 組出 ≤10 列 → results/<帳號>_mva-tiered-v<規則版本>.csv
  make_submission.py check --csv <檔案> [--vcf <work>/norm.vcf.gz] [--scorer evaluation.py]
      格式檢查 → 每個座標回查標準化 VCF → 內建評分模擬 → 官方 evaluation.py 交叉檢查（有下載的話）

內建評分是依官方說明重寫的版本，只用來模擬與交叉檢查；兩者不一致時以官方程式為準。
"""
import argparse
import csv
import importlib.util
import inspect
import os
import re
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import (VcfLookup, default_results, default_rules, default_work, load_rules,  # noqa: E402
                      parse_key, strip_chr, variant_key)

HEADER = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
          "epcr", "finding_type", "notes"]
RANK_POINTS = [(1, 100), (3, 50), (5, 25), (10, 10)]
CHROM_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$")
ALLELE_RE = re.compile(r"^[ACGTN]+$")


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- draft
def route_b(work):
    pairs = read_tsv(os.path.join(work, "trackB_pairs.tsv"))
    info = {}
    for r in read_tsv(os.path.join(work, "trackB_variants.tsv")):
        info[(r["key"], r["gene"])] = r
    out = []
    for p in pairs:
        keys = [k for k in (p["v1"], p["v2"]) if k]
        orig = [info.get((k, p["gene"]), {}).get("ID", "") for k in keys]
        hg = [h for h in (p.get("hgvs1", ""), p.get("hgvs2", "")) if h]
        ev = " | ".join(e for e in (p.get("evidence1", ""), p.get("evidence2", "")) if e)
        note = (f"{p['gene']} (Tier {p['tier']}); route B #{p['rank']} {p['kind']}; phase={p['phase']}; "
                f"{' + '.join(hg)}; evidence: {ev}; orig: {' / '.join(o for o in orig if o)}")
        out.append(dict(src="B", gene=p["gene"], tier=p["tier"], kind=p["kind"], keys=keys, note=note))
    return out


def find_col(header, *names):
    norm = {h.lstrip("#").strip().upper(): h for h in header}
    for n in names:
        if n.upper() in norm:
            return norm[n.upper()]
    return None


def route_a(exdir, proband, lookup):
    """Exomiser 隱性模式（AR）前幾名基因的貢獻變異 → 配對；座標換回標準化 VCF 的寫法。"""
    out = []
    for run in ("exome", "introns"):
        path = os.path.join(exdir, f"{proband}_{run}.variants.tsv")
        if not os.path.exists(path):
            continue
        rows = read_tsv(path)
        if not rows:
            continue
        h = list(rows[0].keys())
        c = {k: find_col(h, *v) for k, v in {
            "rank": ("RANK",), "gene": ("GENE_SYMBOL",), "moi": ("MOI",), "contrib": ("CONTRIBUTING_VARIANT",),
            "contig": ("CONTIG", "CHROM"), "start": ("START", "POS"), "ref": ("REF",), "alt": ("ALT",),
            "gt": ("GENOTYPE", "GT")}.items()}
        if any(v is None for v in c.values()):
            print(f"  ! {os.path.basename(path)} 欄位和預期不同，略過路線 A：缺 {[k for k, v in c.items() if v is None]}")
            continue
        groups = OrderedDict()
        for r in rows:
            if str(r[c["contrib"]]).strip().lower() in ("1", "true", "yes", "y", "t") and r[c["moi"]] in ("AR", "AUTOSOMAL_RECESSIVE"):
                groups.setdefault((int(float(r[c["rank"]])), r[c["gene"]]), []).append(r)
        for (rank, gene), vs in sorted(groups.items()):
            keys = []
            for r in vs:
                hit = lookup.resolve(r[c["contig"]], int(float(r[c["start"]])), r[c["ref"]], r[c["alt"]]) if lookup else None
                if hit:
                    keys.append(variant_key(*hit))
                else:
                    print(f"  ! Exomiser {gene} 的一個變異對不回標準化 VCF，略過")
            hom = len(vs) == 1 and vs[0][c["gt"]].replace("|", "/") == "1/1"
            if len(keys) == 2 or (len(keys) == 1 and hom):
                out.append(dict(src="A", gene=gene, tier="", kind="pair" if len(keys) == 2 else "hom", keys=keys,
                                note=f"{gene}; route A Exomiser {run} #{rank} AR"))
    return out


def hf_user():
    if os.environ.get("HF_USER"):
        return os.environ["HF_USER"]
    try:
        from huggingface_hub import whoami
        return whoami()["name"]
    except Exception:
        return "track1"


def to_cells(keys):
    cells = []
    for k in keys[:2]:
        c, p, r, a = parse_key(k)
        cells += ["chr" + c, str(p), r, a]
    return cells + [""] * (8 - len(cells))


def cmd_draft(a):
    R = load_rules(a.rules)
    lad_p, lad_s = R["submission"]["epcr_primary"], R["submission"]["epcr_secondary"]
    if len(set(lad_p + lad_s)) != len(lad_p + lad_s) or sorted(lad_p + lad_s, reverse=True) != lad_p + lad_s:
        raise SystemExit("ranking_rules.yml 的 EPCR 階梯要嚴格遞減、不能重複")
    B = route_b(a.work)
    vcf = os.path.join(a.work, "norm.vcf.gz")
    lookup = VcfLookup(vcf) if os.path.exists(vcf) else None
    A = route_a(os.path.join(a.work, "exomiser"), a.proband, lookup)

    rows, seen = [], set()

    def add(c, limit):
        fs = frozenset(c["keys"])
        if len(rows) < limit and c["keys"] and fs not in seen:
            seen.add(fs)
            rows.append(c)

    a_pairs = [c for c in B if c["tier"] == "A" and c["kind"] in ("pair", "hom")]
    if a_pairs:
        add(a_pairs[0], 1)                                                    # rank1
        for c in a_pairs[1:]:
            if c["gene"] == a_pairs[0]["gene"]:
                add(c, 4)                                                     # rank2_4：同基因換第二個等位基因
        for c in a_pairs:
            add(c, 4)
    else:
        print("  ! 路線 B 沒有任何 Tier A 配對：第 1 列改用路線 A / 其他候選，務必人工檢查")
    for c in A:
        add(c, 8)                                                             # rank5_8：Exomiser
    for c in B:
        add(c, 8)                                                             # 不足時用路線 B 剩下的
    primary = rows[:]
    secondary = []
    sec_path = os.path.join(a.results, "secondary_findings.tsv")
    if os.path.exists(sec_path):
        for r in read_tsv(sec_path)[:len(lad_s)]:
            k = variant_key(r["chrom"], r["pos"], r["ref"], r["alt"])
            secondary.append(dict(src="SF", gene=r.get("gene", ""), keys=[k],
                                  note=f"secondary finding: {r.get('gene', '')}; {r.get('reason', '')}"))
    if not primary:
        raise SystemExit("沒有任何候選可以提交")

    version = R["version"]
    name = a.name or hf_user()
    out = a.out or os.path.join(a.results, f"{name}_{R['submission']['method_name']}-v{version}.csv")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        w.writerow(HEADER)
        for c, e in list(zip(primary, lad_p)) + list(zip(secondary, lad_s)):
            ft = "secondary" if c["src"] == "SF" else "primary"
            w.writerow([a.proband] + to_cells(c["keys"]) + [f"{e:.2f}", ft, c["note"][:500]])
    print(f"提交草稿：{out}")
    print(f"CSV={out}")
    for i, (c, e) in enumerate(list(zip(primary, lad_p)) + list(zip(secondary, lad_s)), 1):
        print(f"  {i:>2}. EPCR {e:.2f}  {c['src']}  {c['gene']:<10} {' + '.join(c['keys'])}")
    if len(primary[0]["keys"]) != 2:
        print("  ! 第 1 列不是一對變異：答案是一對複合異型合子，第 1 列只放一個變異最多拿一半分數")
    print("下一步：python make_submission.py check --csv " + out)


# ---------------------------------------------------------------- check
def parse_csv(path, proband):
    errs, warns, rows = [], [], []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        data = list(csv.reader(fh))
    if not data or [h.strip() for h in data[0]] != HEADER:
        errs.append(f"表頭必須是：{','.join(HEADER)}")
        return rows, errs, warns
    body = [r for r in data[1:] if any(x.strip() for x in r)]
    if not 1 <= len(body) <= 10:
        errs.append(f"要有 1–10 列候選，目前 {len(body)} 列")
    for i, r in enumerate(body, 1):
        if len(r) != len(HEADER):
            errs.append(f"第 {i} 列有 {len(r)} 欄，要 {len(HEADER)} 欄（notes 有逗號要用雙引號）")
            continue
        d = dict(zip(HEADER, [x.strip() for x in r]))
        if d["proband_id"] != proband:
            errs.append(f"第 {i} 列 proband_id 是 {d['proband_id']!r}，要是 {proband}")
        vs = []
        for n in ("1", "2"):
            cell = [d[f"chrom_{n}"], d[f"pos_{n}"], d[f"ref_{n}"], d[f"alt_{n}"]]
            if not any(cell):
                if n == "1":
                    errs.append(f"第 {i} 列沒有第 1 個變異")
                continue
            if not all(cell):
                errs.append(f"第 {i} 列第 {n} 個變異欄位不完整：{cell}")
                continue
            ch, pos, ref, alt = cell
            if not CHROM_RE.match(ch):
                errs.append(f"第 {i} 列 chrom_{n} = {ch!r}：要寫成 chr15 這種格式（VCF 裡沒有 chr，提交要加）")
            if not pos.isdigit() or int(pos) < 1:
                errs.append(f"第 {i} 列 pos_{n} 不是正整數：{pos!r}")
                continue
            if not ALLELE_RE.match(ref) or not ALLELE_RE.match(alt) or ref == alt:
                errs.append(f"第 {i} 列第 {n} 個變異 REF/ALT 不合法（大寫 ACGTN、REF≠ALT）：{ref}>{alt}")
            vs.append((ch, int(pos), ref, alt))
        try:
            e = float(d["epcr"])
            if not 0 < e <= 1:
                errs.append(f"第 {i} 列 EPCR 要在 (0, 1]：{e}")
        except ValueError:
            errs.append(f"第 {i} 列 EPCR 不是數字：{d['epcr']!r}")
            e = float("nan")
        if d["finding_type"] not in ("primary", "secondary"):
            errs.append(f"第 {i} 列 finding_type 要是 primary 或 secondary：{d['finding_type']!r}")
        rows.append(dict(i=i, variants=frozenset(vs), epcr=e, ft=d["finding_type"], n=len(vs)))
    eps = [r["epcr"] for r in rows]
    if len(set(eps)) != len(eps):
        errs.append("EPCR 有重複：每一列都要不同，尤其不能和第 1 名同分（會拉低 F-max）")
    if rows:
        top = max(rows, key=lambda r: (r["epcr"], -r["i"]))
        if top["n"] != 2:
            warns.append("EPCR 最高的那一列不是一對變異：答案是一對複合異型合子，單一變異最多拿一半分數")
        if top["ft"] != "primary":
            warns.append("EPCR 最高的那一列不是 primary")
    return rows, errs, warns


def rank_points(rank):
    for lim, pts in RANK_POINTS:
        if rank <= lim:
            return pts
    return 0


def simulate(rows, truth):
    """依官方說明：依 EPCR 由高到低排名（同分照原順序）；完全命中 = 同一列剛好是答案的那一組變異；
    只命中其中一個變異 = 一半分數。F-max：每個不同的 EPCR 門檻下，比對 ≥ 門檻的所有變異。"""
    ranked = sorted(rows, key=lambda r: (-r["epcr"], r["i"]))
    best = (0, "none", None)
    for k, r in enumerate(ranked, 1):
        if r["variants"] == truth:
            cand = (rank_points(k), "full", k)
        elif r["variants"] & truth:
            cand = (rank_points(k) / 2, "partial", k)
        else:
            continue
        if cand[0] > best[0]:
            best = cand
    fmax = 0.0
    for t in sorted({r["epcr"] for r in rows}, reverse=True):
        pred = set().union(*[r["variants"] for r in rows if r["epcr"] >= t])
        tp = len(pred & truth)
        if tp:
            p, rc = tp / len(pred), tp / len(truth)
            fmax = max(fmax, 2 * p * rc / (p + rc))
    return best[0], round(fmax, 4), best[1], best[2]


def official(path, csv_path, proband, truth_row_index):
    spec = importlib.util.spec_from_file_location("mva_official_evaluation", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sub = mod.load_submission(csv_path)
    prow = sub[proband] if isinstance(sub, dict) else sub
    srt = sorted(prow, key=lambda r: (-float(getattr(r, "epcr", 0)), getattr(r, "rank", 0)))
    truth = srt[truth_row_index].variants
    n = len(inspect.signature(mod.score_proband).parameters)
    res = mod.score_proband(proband, prow, truth) if n >= 3 else mod.score_proband(prow, truth)
    get = (lambda k: res.get(k)) if isinstance(res, dict) else (lambda k: getattr(res, k, None))
    pts = next((get(k) for k in ("rank_points", "points", "score") if get(k) is not None), None)
    fm = next((get(k) for k in ("f_max", "fmax", "F_max") if get(k) is not None), None)
    return pts, fm, res


def cmd_check(a):
    rows, errs, warns = parse_csv(a.csv, a.proband)
    print(f"檢查 {a.csv}：{len(rows)} 列")
    if rows and a.vcf and os.path.exists(a.vcf):
        lk = VcfLookup(a.vcf)
        miss = [(r["i"], v) for r in rows for v in sorted(r["variants"]) if not lk.exact(strip_chr(v[0]), v[1], v[2], v[3])]
        for i, v in miss:
            errs.append(f"第 {i} 列的 {v[0]}:{v[1]} {v[2]}>{v[3]} 在標準化 VCF 裡找不到（抄錯，或寫法和 norm.vcf.gz 不同）")
        if not miss:
            print("  ✓ 每個變異都在標準化 VCF 裡找得到，REF/ALT 完全一致")
    else:
        warns.append("沒有回查 VCF（找不到 norm.vcf.gz）")
    for w in warns:
        print("  ! " + w)
    for e in errs:
        print("  ✗ " + e)
    if errs:
        print("格式檢查沒過，先修正再提交。")
        sys.exit(1)
    print("  ✓ 格式檢查通過")

    print("\n本機模擬（假設某一列就是答案時會拿幾分）")
    ranked = sorted(rows, key=lambda r: (-r["epcr"], r["i"]))
    sims = []
    for k, r in enumerate(ranked, 1):
        if r["n"] == 2:
            pts, fm, kind, at = simulate(rows, r["variants"])
            sims.append((k, pts, fm))
            print(f"  答案 = 第 {k} 名（EPCR {r['epcr']:.2f}）→ Rank points {pts:g}、F-max {fm:.3f}（{kind} @ {at}）")
    if ranked and ranked[0]["n"] == 2:
        one = sorted(ranked[0]["variants"])[0]
        pts, fm, kind, at = simulate(rows, frozenset({one, ("chr0", 1, "A", "C")}))
        print(f"  第 1 名只對一個變異（另一個不在清單裡）→ Rank points {pts:g}、F-max {fm:.3f}（{kind} @ {at}）")

    if a.scorer and os.path.exists(a.scorer):
        print(f"\n官方 evaluation.py 交叉檢查（{a.scorer}）")
        try:
            for k, pts, fm in sims[:2]:
                o_pts, o_fm, raw = official(a.scorer, a.csv, a.proband, k - 1)
                same = o_pts is not None and abs(float(o_pts) - pts) < 1e-6 and (o_fm is None or abs(float(o_fm) - fm) < 1e-3)
                print(f"  答案 = 第 {k} 名：官方 Rank points {o_pts}、F-max {o_fm} → {'和本機一致' if same else '和本機不同，以官方為準'}")
                if o_pts is None:
                    print(f"    官方回傳：{raw!r}")
            print("  ✓ 官方 load_submission 讀得進這個檔")
        except Exception as e:  # 官方程式介面改了也不要讓整個檢查失敗
            print(f"  ! 官方評分程式執行失敗：{type(e).__name__}: {e}")
            print("    先確認 CSV 格式；若官方程式介面改了，打開 evaluation.py 看 load_submission / score_proband 的用法")
    else:
        print("\n（沒有官方 evaluation.py：跑 02_fetch_resources.sh scorer 下載後再檢查一次）")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("draft")
    d.add_argument("--work", default=default_work())
    d.add_argument("--results", default=default_results())
    d.add_argument("--rules", default=default_rules())
    d.add_argument("--name", help="HF 帳號（檔名用），不給就自動讀 hf 登入的帳號")
    d.add_argument("--out")
    d.add_argument("--proband", default=os.environ.get("PROBAND", "PROBAND01"))
    c = sub.add_parser("check")
    c.add_argument("--csv", required=True)
    c.add_argument("--vcf", default=os.path.join(default_work(), "norm.vcf.gz"))
    c.add_argument("--scorer", default=os.path.expanduser(os.environ.get("SCORER") or "~/mva-hackathon/resources/scorer/evaluation.py"))
    c.add_argument("--proband", default=os.environ.get("PROBAND", "PROBAND01"))
    a = ap.parse_args()
    cmd_draft(a) if a.cmd == "draft" else cmd_check(a)


if __name__ == "__main__":
    main()
