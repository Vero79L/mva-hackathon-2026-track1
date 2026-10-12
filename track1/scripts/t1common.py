"""Track 1 共用工具：路徑、規則、基因分層、變異寫法、VCF 查詢、HPO。

只用 Python 標準函式庫 + PyYAML（pandas 只在 rank_panel.py 用）。
"""
import gzip
import os
import re
import shutil
import subprocess

import yaml

T1_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONF_DIR = os.path.join(T1_DIR, "config")


def env_path(name, default):
    """環境變數（由 lib.sh 從 config/paths.env 匯出）優先，否則用預設值。"""
    return os.path.expanduser(os.environ.get(name) or default)


def default_work():
    return env_path("WORK", "~/mva-hackathon/work/track1")


def default_results():
    return env_path("RESULTS", "~/mva-hackathon/results/track1")


def default_rules():
    return os.path.join(CONF_DIR, "ranking_rules.yml")


# ---------- 規則 ----------
POINT_KEYS = ["PVS1", "clinvar_plp_2star", "clinvar_plp_1star", "known_allele", "spliceai_strong",
              "spliceai_supporting", "pp3_sift_polyphen", "tierA_missense_rare", "tierA_noncoding_rare",
              "pm2_supporting", "pp4_tierA"]
TIEBREAKS = {"weaker_impact", "phase_trans", "sum_points", "quality", "rarer"}


def load_rules(path):
    with open(path, encoding="utf-8") as fh:
        r = yaml.safe_load(fh)
    need = ["version", "genes", "frequency", "spliceai", "points", "pairing", "merge", "submission"]
    missing = [k for k in need if k not in r]
    if missing:
        raise SystemExit(f"ranking_rules.yml 缺少段落：{', '.join(missing)}")
    missing = [k for k in POINT_KEYS if k not in r["points"]]
    if missing:
        raise SystemExit(f"ranking_rules.yml 的 points 缺少：{', '.join(missing)}")
    bad = [t for t in r["pairing"].get("tiebreak", []) if t not in TIEBREAKS]
    if bad:
        raise SystemExit(f"ranking_rules.yml 的 pairing.tiebreak 有不認識的項目：{bad}（可用：{sorted(TIEBREAKS)}）")
    r.setdefault("known_alleles", [])
    r.setdefault("quality", {})
    r["quality"].setdefault("min_gq", 20)
    r["quality"].setdefault("min_dp", 8)
    r["frequency"].setdefault("pm2_max_af", 1e-4)
    r["pairing"].setdefault("single_het_factor", 0.5)
    r["pairing"].setdefault("exclude_cis", True)
    for ka in r["known_alleles"]:
        parse_key(ka["variant"])           # 寫法錯就在這裡報錯
    return r


def gene_sets(rules):
    g = rules["genes"]
    return set(g["tier_a"]), set(g["tier_b"]), dict(g.get("aliases") or {})


# ---------- 變異寫法 ----------
def strip_chr(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def variant_key(chrom, pos, ref, alt):
    return f"{strip_chr(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def parse_key(k):
    parts = str(k).split(":")
    if len(parts) != 4 or not parts[1].isdigit() or not re.fullmatch(r"[ACGTN]+", parts[2].upper()) \
            or not re.fullmatch(r"[ACGTN]+", parts[3].upper()):
        raise ValueError(f"變異寫法要是 CHROM:POS:REF:ALT（例如 15:40117088:G:A）：{k}")
    return strip_chr(parts[0]), int(parts[1]), parts[2].upper(), parts[3].upper()


def minrep(pos, ref, alt):
    """最精簡寫法：去掉共同的尾端，再去掉共同的開頭（用來比對不同工具的 indel 寫法）。"""
    ref = "" if ref in ("-", ".", None) else str(ref).upper()
    alt = "" if alt in ("-", ".", None) else str(alt).upper()
    pos = int(pos)
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


# ---------- VCF 查詢（有 bcftools + 索引就用 bcftools；否則用 Python 讀，較慢） ----------
class VcfLookup:
    def __init__(self, vcf):
        self.vcf = vcf
        self.bcftools = shutil.which("bcftools")
        self.indexed = os.path.exists(vcf + ".tbi") or os.path.exists(vcf + ".csi")
        self._mem = None

    def _open(self):
        return gzip.open(self.vcf, "rt") if self.vcf.endswith(".gz") else open(self.vcf)

    def records(self, chrom, start, end):
        """回傳與區間重疊的紀錄 [(chrom, pos, ref, alt, id)]，chrom 不帶 chr。"""
        chrom = strip_chr(chrom)
        if self.bcftools and self.indexed:
            out = subprocess.run([self.bcftools, "query", "-r", f"{chrom}:{max(1, start)}-{end}",
                                  "-f", r"%CHROM\t%POS\t%REF\t%ALT\t%ID\n", self.vcf],
                                 capture_output=True, text=True, check=True).stdout
            rows = [ln.split("\t") for ln in out.splitlines() if ln]
            return [(strip_chr(c), int(p), r, a, i) for c, p, r, a, i in rows]
        if self._mem is None and os.path.getsize(self.vcf) < 50_000_000:
            self._mem = list(self._iter())
        it = self._mem if self._mem is not None else self._iter()
        return [x for x in it if x[0] == chrom and x[1] <= end and x[1] + len(x[2]) - 1 >= start]

    def _iter(self):
        with self._open() as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                f = ln.split("\t", 5)
                for alt in f[4].split(","):
                    yield strip_chr(f[0]), int(f[1]), f[3].upper(), alt.upper(), f[2]

    def exact(self, chrom, pos, ref, alt):
        return any(p == int(pos) and r == ref.upper() and a == alt.upper()
                   for _, p, r, a, _ in self.records(chrom, int(pos), int(pos)))

    def resolve(self, chrom, pos, ref, alt, window=10):
        """把其他工具的寫法對回這個 VCF 的標準化寫法；找不到回傳 None。"""
        target = minrep(pos, ref, alt)
        for c, p, r, a, _ in self.records(chrom, int(pos) - window, int(pos) + window):
            if minrep(p, r, a) == target:
                return c, p, r, a
        return None


def vcf_samples(vcf):
    if shutil.which("bcftools"):
        out = subprocess.run(["bcftools", "query", "-l", vcf], capture_output=True, text=True, check=True).stdout
        return [s for s in out.split() if s]
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#CHROM"):
                return ln.rstrip("\n").split("\t")[9:]
    return []


# ---------- HPO ----------
HPO_RE = re.compile(r"\bHP[:_](\d{7})\b")


def parse_obo(path):
    """hp.obo → {'terms': {id: {name, obsolete, replaced_by}}, 'alt': {alt_id: id}, 'version': str}"""
    terms, alt, version, cur = {}, {}, "", None
    if not path or not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith("data-version:") and not version:
                version = ln.split(":", 1)[1].strip()
            elif ln == "[Term]":
                cur = {"id": None, "name": "", "obsolete": False, "replaced_by": None}
            elif ln.startswith("[") and ln.endswith("]"):
                cur = None
            elif cur is not None and ":" in ln:
                k, v = ln.split(":", 1)
                v = v.strip()
                if k == "id":
                    cur["id"] = v
                    terms[v] = cur
                elif k == "name":
                    cur["name"] = v
                elif k == "is_obsolete" and v == "true":
                    cur["obsolete"] = True
                elif k == "replaced_by":
                    cur["replaced_by"] = v
                elif k == "alt_id" and cur["id"]:
                    alt[v] = cur["id"]
    return {"terms": terms, "alt": alt, "version": version}


def resolve_hpo(hid, obo):
    """回傳 (最終 ID, 名稱, 說明)。最終 ID 為 None 代表應該丟掉。"""
    if obo is None:
        return hid, "", "沒有 hp.obo，未檢查"
    t = obo["terms"].get(hid)
    if t is None and hid in obo["alt"]:
        new = obo["alt"][hid]
        return new, obo["terms"][new]["name"], f"{hid} 是 {new} 的舊編號（alt_id），已換成新編號"
    if t is None:
        return hid, "", f"{hid} 不在 hp.obo 裡（打錯？），保留但請人工確認"
    if t["obsolete"]:
        rb = t.get("replaced_by")
        if rb and rb in obo["terms"]:
            return rb, obo["terms"][rb]["name"], f"{hid} 已過時，換成 {rb}"
        return None, t["name"], f"{hid} 已過時且沒有替代編號，已移除"
    return hid, t["name"], ""
