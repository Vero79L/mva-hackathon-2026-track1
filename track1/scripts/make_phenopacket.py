#!/usr/bin/env python3
"""從 HPO 清單產生 Exomiser 用的 phenopacket（v2 YAML）。

  make_phenopacket.py --hpo <results>/hpo_terms.tsv --vcf <work>/norm.vcf.gz --sex MALE|FEMALE|UNKNOWN \
                      --out <work>/exomiser/proband.phenopacket.yml [--obo hp.obo] [--id PROBAND01]

subject.id 自動取 VCF 裡的樣本名稱（Exomiser 要求兩者一致）。輸出含表現型，放 work/，不進 repo。
"""
import argparse
import datetime
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t1common import HPO_RE, parse_obo, resolve_hpo, vcf_samples  # noqa: E402

SEX = {"MALE": "MALE", "M": "MALE", "FEMALE": "FEMALE", "F": "FEMALE", "UNKNOWN": "UNKNOWN_SEX", "": "UNKNOWN_SEX"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hpo", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sex", default="UNKNOWN")
    ap.add_argument("--out", required=True)
    ap.add_argument("--obo")
    ap.add_argument("--id", default=os.environ.get("PROBAND", "PROBAND01"))
    ap.add_argument("--sample", help="不給就讀 VCF 的樣本名稱")
    a = ap.parse_args()

    sex = SEX.get(a.sex.strip().upper())
    if sex is None:
        raise SystemExit(f"--sex 要是 MALE、FEMALE 或 UNKNOWN：{a.sex}")
    if sex == "UNKNOWN_SEX":
        print("  ! 性別是 UNKNOWN：Exomiser 仍可執行，但 X 連鎖遺傳模式的判斷會變弱。確認後在 config/local.env 設 SEX=MALE 或 FEMALE。")

    samples = [a.sample] if a.sample else vcf_samples(a.vcf)
    if len(samples) != 1:
        raise SystemExit(f"VCF 裡有 {len(samples)} 個樣本，用 --sample 指定：{samples}")

    obo = parse_obo(a.obo)
    feats, seen = [], set()
    with open(a.hpo, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            m = HPO_RE.search(f[0]) if f and f[0] else None
            if not m:
                continue
            hid, label = "HP:" + m.group(1), (f[1].strip() if len(f) > 1 else "")
            new, name, note = resolve_hpo(hid, obo)
            if note:
                print("  ! " + note)
            if new is None or new in seen:
                continue
            seen.add(new)
            feats.append({"type": {"id": new, "label": name or label or new}})
    if not feats:
        raise SystemExit(f"{a.hpo} 裡沒有有效的 HPO ID")

    pp = {
        "id": a.id,
        "subject": {"id": samples[0], "sex": sex},
        "phenotypicFeatures": feats,
        "metaData": {
            "created": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "createdBy": "mva-track1",
            "resources": [{"id": "hp", "name": "human phenotype ontology",
                           "url": "http://purl.obolibrary.org/obo/hp.owl",
                           "version": (obo or {}).get("version") or "unknown",
                           "namespacePrefix": "HP", "iriPrefix": "http://purl.obolibrary.org/obo/HP_"}],
            "phenopacketSchemaVersion": "2.0",
        },
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        yaml.safe_dump(pp, fh, allow_unicode=True, sort_keys=False)
    print(f"phenopacket：{len(feats)} 個 HPO terms，subject.id = {samples[0]}，sex = {sex} → {a.out}")


if __name__ == "__main__":
    main()
