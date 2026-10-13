"""測試用的「官方評分程式介面」替身——不是官方檔案。

只模仿 make_submission.py 會呼叫的兩個函式（load_submission、score_proband）和回傳物件的屬性，
用來測試介面轉接的程式碼。真正提交前，08_submission.sh 會用 02 下載的官方 evaluation.py。

from __future__ import annotations 要保留：官方檔也有這行，它讓 @dataclass 的欄位註解變成字串、
走進 dataclasses._is_type 查 sys.modules 的路——official() 忘了登記模組時就會在這裡 AttributeError。
"""
from __future__ import annotations

import csv
from dataclasses import dataclass


@dataclass
class Row:
    rank: int
    epcr: float
    finding_type: str
    variants: frozenset


@dataclass
class Result:
    rank_points: float
    f_max: float


def load_submission(path):
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    parsed = []
    for i, r in enumerate(rows):
        vs = set()
        for n in ("1", "2"):
            if r[f"chrom_{n}"]:
                vs.add((r[f"chrom_{n}"], int(r[f"pos_{n}"]), r[f"ref_{n}"], r[f"alt_{n}"]))
        parsed.append((i, r["proband_id"], float(r["epcr"]), r["finding_type"], frozenset(vs)))
    for pid in {p[1] for p in parsed}:
        mine = sorted([p for p in parsed if p[1] == pid], key=lambda p: (-p[2], p[0]))
        out[pid] = [Row(k, p[2], p[3], p[4]) for k, p in enumerate(mine, 1)]
    return out


def score_proband(proband_id, rows, truth):
    truth = frozenset(truth)
    pts = 0
    for r in rows:
        base = 100 if r.rank == 1 else 50 if r.rank <= 3 else 25 if r.rank <= 5 else 10 if r.rank <= 10 else 0
        if r.variants == truth:
            pts = max(pts, base)
        elif r.variants & truth:
            pts = max(pts, base / 2)
    fmax = 0.0
    for t in sorted({r.epcr for r in rows}, reverse=True):
        pred = set().union(*[r.variants for r in rows if r.epcr >= t])
        tp = len(pred & truth)
        if tp:
            p, rc = tp / len(pred), tp / len(truth)
            fmax = max(fmax, 2 * p * rc / (p + rc))
    return Result(pts, fmax)
