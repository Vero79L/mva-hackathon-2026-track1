# Track 1 Methods — Vero7878

## 1. Summary

Top-ranked candidate: **<填: 基因>** compound heterozygous **<填: HGVS c./p. allele 1>** / **<填: HGVS allele 2>**
(GRCh38 <填: chr:pos:ref:alt> and <填: chr:pos:ref:alt>), submitted at EPCR 0.95.
<填: 一句話理由，例如 both routes rank the same gene first / a truncating allele paired with a damaging missense allele.>

## 2. Data

- Organizer-provided WGS VCF (`WGS_*.vcf.gz`; Sentieon 202308.02 + GATK VariantFiltration; GRCh38 no-alt analysis set
  with hs38d1 decoys, no `chr` prefix); single sample PROBAND01.
- Organizer-provided HPO terms (<填: n> terms; obsolete IDs replaced using hp.obo release 2026-09-01).
- No FASTQ/BAM re-analysis in Phase 1. No parental data available.

## 3. Pipeline and software

| Step | Tool / version | Notes |
|---|---|---|
| QC | bcftools 1.24 | Genome-wide aggregates in `qc/qc_summary.txt`; sex inferred from chrX non-PAR het/hom-alt ratio and chrY calls |
| Normalization | bcftools norm (1.24) | Multi-allelic sites split, indels left-aligned against the matching reference, `*` alleles removed, original representation kept in ID |
| Annotation | Ensembl VEP 116.2 (offline cache 116; gnomAD v4.1, MANE v1.5), ClinVar 2026-10-04 | 13-gene panel ±10 kb plus known distal alleles; `--everything --hgvsg` |
| Splicing | SpliceAI 1.3.1 (TensorFlow 2.15), ±500 bp | Max of the four delta scores, counted only for the matching gene |
| Route A | Exomiser 15.1.0, data 2602 (Java 21.0.12.1) | `exome` preset and exome-with-introns preset |
| Route B | `scripts/rank_panel.py` + `config/ranking_rules.yml` v1 | See section 4 |

Python 3.11.17, pandas 3.0.6; macOS 26.5 (arm64). Exact environments: `env/*.yml`, `env/tools.txt`, `env/java.txt`,
`env/exomiser.txt`, `env/clinvar_version.txt`.

## 4. Pre-registered rules (route B)

Rules were committed at **17ea747** on **2026-10-12 15:39:16 +0800**, before the first ranking run on
**2026-10-13 16:20:41** (`docs/preregistration/first_rank_run.txt`, `docs/preregistration/run_manifest.tsv`;
each run records the rules commit, version and SHA-256).

- Gene tiers: Tier A = BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5 (MVA1–6); Tier B = BUB1, BUB3, MAD1L1, MAD2L1BP,
  CEP192, KNL1, ZWINT.
- Allele points (simplified ACMG/AMP Bayesian scale): PVS1 8 (Tier A LoF); ClinVar P/LP ≥2★ 4, 1★ 2; known
  hypomorphic allele 4 (max with ClinVar, not additive); SpliceAI ≥0.5 2, ≥0.2 1; PP3 (SIFT and PolyPhen damaging) 1;
  rare Tier A missense/in-frame 1; rare Tier A non-coding 1; PM2_supporting (MAX_AF < 1e-4) 1; PP4 (Tier A) 1.
- Qualifying allele: MAX_AF ≤ 0.01 (or ClinVar P/LP, or a listed known allele) + at least one functional evidence item.
- Pairing: pair score = points of the weaker allele; pairs proven *cis* by PGT/PID are excluded; ties broken by
  weaker-allele VEP impact, phase-confirmed *trans*, total points, call quality, then lower allele frequency.
- Merge into ≤10 rows: rank 1 = best Tier A pair; ranks 2–4 = alternative pairs in the same gene; ranks 5–8 =
  Exomiser recessive-model pairs; ranks 9–10 = reviewed secondary findings. EPCR 0.95, 0.60, 0.45, 0.35, 0.25, 0.18,
  0.12, 0.08 (secondary 0.05, 0.04).

## 5. Results

- Route A: <填: top genes and the rank of each MVA gene, from results/track1/exomiser_summary.md>.
- Route B: <填: top pairs; evidence for each allele; phase (trans / unphased → presumed trans)>.
- Agreement between routes: <填: yes/no>.
- Submitted CSV: 8 rows; format check, VCF back-check and self-scoring (local simulation cross-checked against the
  official `evaluation.py`) all passed before submission.
- Leaderboard feedback for submission #1 (2026-10-14): Rank points **100.0 / 100**, F-max **1.000**
  (EPCR threshold 0.95), **full match at rank 1**.

## 6. Limitations

- No parental data: unphased pairs are presumed *trans*.
- Short reads miss structural variants, mobile-element (Alu) insertions and some deep-intronic events; SV/CNV calls
  overlapping the panel were reviewed (no SV/CNV VCF was provided in the dataset).
- The ACMG/AMP-like points are a heuristic, not a clinical classification.
- Route B is restricted to the 13-gene panel; genome-wide ranking relies on Exomiser.

## 7. LLM / AI usage

- **Code and documentation**: Claude (Anthropic), used on a Max subscription (consumer terms) through claude.ai and
  Claude Code, wrote the pipeline code, tests and documentation. Claude Code ran with user-level read-deny rules and
  an OS-level sandbox covering `gated-data/`, `work/`, `results/` and `track1/qc/`; it never read the raw data or any
  file in those folders.
- **Running the pipeline**: every step that reads patient data (03–08) was started by the participant in a separate
  terminal via `scripts/supervise.sh`, a plain bash script with no LLM component. No Agent SDK, API-authenticated
  model, or any other model was used to run or diagnose the pipeline in Phase 1.
- **Handoff**: the designed channel from the pipeline back to the assistant is `handoff/STATUS.md` — step names,
  pass/fail, timings and an error class, written by code with programmatic masking (coordinates, digit runs, alleles,
  HGVS, rsIDs, HPO IDs, genotypes, sample names, panel gene names, sex, data file names).
- **Model exposure to candidate-level results**: after the pipeline had completed, the participant pasted terminal
  output of the finished runs into the Claude Code conversation. This exposed candidate-level results — the ranked
  candidate pairs (coordinates and HGVS), gene-level Exomiser ranks, the VCF sample ID and the inferred sex — to the
  consumer-terms assistant. The assistant did not access the data folders at any point, and the analysis code and
  ranking rules were committed before the first ranking run (`docs/preregistration/`).
- **Data processing settings**: claude.ai "Help improve Claude" (model training) is **off** on the account used; the
  same setting applies to Claude Code, which signs in with the same account.

## 8. Acknowledgement

<填: Official Rules 頁面最下方要求的 acknowledgement 原文>
