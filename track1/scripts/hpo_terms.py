#!/usr/bin/env python3
"""從比賽資料夾的表現型檔抽出 HPO ID，對照 hp.obo 補上名稱、換掉過時編號。

  hpo_terms.py extract --data <gated-data> --out <results>/hpo_terms.tsv [--obo hp.obo]
  hpo_terms.py check   --tsv <hpo_terms.tsv> [--obo hp.obo]     檢查手動整理的清單

輸出 TSV 三欄：hpo_id、label、note。這是孩子的表現型，放 results/，不進 repo。
找不到任何 HPO ID 時（例如表現型只有文字描述），會建一個空白範本並以代碼 2 結束。
"""
import argparse
import os
import re
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import HPO_RE, parse_obo, resolve_hpo  # noqa: E402

TEXT_EXT = {".txt", ".tsv", ".csv", ".json", ".md", ".yml", ".yaml", ".xml", ".html", ".htm", ".tab"}
SKIP_EXT = {".gz", ".bgz", ".tbi", ".csi", ".bam", ".bai", ".cram", ".crai", ".fastq", ".fq", ".bcf"}
LABEL_RE = re.compile(r"HP[:_]\d{7}\s*[\]\)]?\s*[:\-–—|,\t ]\s*([A-Za-z][A-Za-z0-9 ,'()/\-]{2,80})")


def read_text(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in TEXT_EXT:
        with open(path, encoding="utf-8", errors="ignore") as fh:
            return fh.read()
    if ext == ".xlsx":
        try:
            import openpyxl
        except ImportError:
            print(f"  （略過 {os.path.basename(path)}：沒有 openpyxl）")
            return ""
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        return "\n".join("\t".join("" if c is None else str(c) for c in row)
                         for ws in wb.worksheets for row in ws.iter_rows(values_only=True))
    if ext == ".docx":
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8", "ignore")
        return re.sub(r"<[^>]+>", " ", xml.replace("</w:p>", "\n"))
    if ext == ".pdf":
        try:
            from pypdf import PdfReader
        except ImportError:
            print(f"  （略過 {os.path.basename(path)}：沒有 pypdf；可手動整理）")
            return ""
        return "\n".join((p.extract_text() or "") for p in PdfReader(path).pages)
    return ""


def extract(data_dir):
    found, sources = {}, {}
    for root, _, files in os.walk(data_dir):
        for fn in sorted(files):
            path = os.path.join(root, fn)
            low = fn.lower()
            if any(low.endswith(e) for e in SKIP_EXT) or ".fastq" in low or os.path.getsize(path) > 50_000_000:
                continue
            try:
                text = read_text(path)
            except Exception as e:  # 單一檔案讀不了不影響其他檔
                print(f"  （讀不了 {fn}：{e}）")
                continue
            for line in text.splitlines():
                for m in HPO_RE.finditer(line):
                    hid = "HP:" + m.group(1)
                    sources.setdefault(hid, set()).add(os.path.relpath(path, data_dir))
                    if not found.get(hid):
                        lm = LABEL_RE.search(line[m.start():])
                        found[hid] = lm.group(1).strip() if lm else found.get(hid, "")
    return found, sources


def write(rows, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("hpo_id\tlabel\tnote\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def clean(terms, obo):
    rows, seen = [], set()
    for hid, label in terms:
        new, name, note = resolve_hpo(hid, obo)
        if note:
            print("  ! " + note)
        if new is None or new in seen:
            continue
        seen.add(new)
        rows.append((new, name or label, note))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--data", required=True)
    e.add_argument("--out", required=True)
    e.add_argument("--obo")
    c = sub.add_parser("check")
    c.add_argument("--tsv", required=True)
    c.add_argument("--obo")
    a = ap.parse_args()
    obo = parse_obo(a.obo)
    if obo is None:
        print("  （沒有 hp.obo：不檢查過時編號。先跑 02_fetch_resources.sh hpo）")

    if a.cmd == "extract":
        found, sources = extract(a.data)
        if not found:
            write([], a.out)
            print(f"在 {a.data} 沒有找到任何 HPO ID。打開表現型檔，把 HPO ID 和名稱手動填進 {a.out}（Tab 分隔），"
                  "再跑 hpo_terms.py check。")
            sys.exit(2)
        rows = clean(list(found.items()), obo)
        write(rows, a.out)
        print(f"{len(rows)} 個 HPO terms → {a.out}")
        for src in sorted({s for v in sources.values() for s in v}):
            print(f"  來源：{src}（{sum(1 for v in sources.values() if src in v)} 個）")
        return

    terms = []
    with open(a.tsv, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            m = HPO_RE.search(f[0]) if f and f[0] else None
            if m:
                terms.append(("HP:" + m.group(1), f[1].strip() if len(f) > 1 else ""))
            elif f and f[0].strip() and not f[0].startswith(("#", "hpo_id")):
                print(f"  ! 看不懂這一行，略過：{ln.strip()}")
    rows = clean(terms, obo)
    print(f"{len(rows)} 個有效 HPO terms")


if __name__ == "__main__":
    main()
