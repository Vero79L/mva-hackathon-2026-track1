# Track 1 Methods — <HF username or team name>

<!--
草稿：第一階段寫短版，第二階段擴寫。尖括號 <…> 是要填的地方。
只能寫你提名的少數變異；不能貼任何全基因體的基因型表。
官方 methods 模板（methods_description_form.xlsx）的必填欄位要逐一對到下面的段落。
-->

## 1. Summary

Top-ranked candidate: **<GENE>** compound heterozygous **<HGVS c./p. allele 1>** / **<HGVS allele 2>**
(GRCh38 <chr:pos:ref:alt> and <chr:pos:ref:alt>), submitted at EPCR 0.95.
<One sentence on why: e.g. both routes agree / truncating allele plus a known hypomorphic allele.>

## 2. Data

- Organizer-provided WGS VCF (`WGS_*.vcf.gz`; Sentieon 202308.02 + GATK VariantFiltration; GRCh38 no-alt analysis set
  with hs38d1 decoys, no `chr` prefix); single sample PROBAND01.
- Organizer-provided HPO terms (<n> terms; obsolete IDs replaced using hp.obo <version>).
- No FASTQ/BAM re-analysis in Phase 1. No parental data available.

## 3. Pipeline and software

| Step | Tool / version | Notes |
|---|---|---|
| QC | bcftools <version> | Genome-wide aggregates in `qc/qc_summary.txt`; sex inferred from chrX non-PAR het/hom-alt ratio and chrY calls |
| Normalization | bcftools norm | Multi-allelic sites split, indels left-aligned against the matching reference, `*` alleles removed, original representation kept in ID |
| Annotation | Ensembl VEP 116 (offline cache; gnomAD v4.1, MANE v1.5), ClinVar <fileDate> | 13-gene panel ±10 kb plus known distal alleles; `--everything --hgvsg` |
| Splicing | SpliceAI 1.3.1 (TensorFlow 2.15), ±500 bp | Max of the four delta scores, counted only for the matching gene |
| Route A | Exomiser 15.1.0, data 2602 | `exome` preset and exome-with-introns preset |
| Route B | `scripts/rank_panel.py` + `config/ranking_rules.yml` v<version> | See section 4 |

Exact environments: `env/*.yml`, `env/tools.txt`, `env/java.txt`, `env/exomiser.txt`, `env/clinvar_version.txt`.

## 4. Pre-registered rules (route B)

Rules were committed at **<commit hash>** on **<date time>**, before the first ranking run on **<date time>**
(`results/track1/first_rank_run.txt`, `run_manifest.tsv`).

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

- Route A: <top genes and the rank of each MVA gene, from results/track1/exomiser_summary.md>.
- Route B: <top pairs; evidence for each allele; phase (trans / unphased → presumed trans)>.
- Agreement between routes: <yes/no>.
- Leaderboard feedback for submission #<n>: Rank points <…>, F-max <…>, <full / partial / no match>.

## 6. Limitations

- No parental data: unphased pairs are presumed *trans*.
- Short reads miss structural variants, mobile-element (Alu) insertions and some deep-intronic events; SV/CNV calls
  overlapping the panel were reviewed manually.
- The ACMG/AMP-like points are a heuristic, not a clinical classification.
- Route B is restricted to the 13-gene panel; genome-wide ranking relies on Exomiser.

## 7. LLM / AI usage

<Confirm each line against what actually happened, then delete this note.>

- **Code and documentation**: Claude (Anthropic) on a Max subscription, through claude.ai and Claude Code. Neither
  had access to patient data: the chat was connected only to the code folders, and Claude Code ran with read-deny
  rules and an OS-level sandbox over the data folders.
- **Running the pipeline**: steps that read patient data ran locally as plain scripts. When a step failed, an agent
  built on the Claude Agent SDK, authenticated with an API key under Anthropic's Commercial Terms (API usage paid from
  subscriber API credits), diagnosed environment problems on the same machine. It could not edit the analysis code or
  the pre-registered rules.
- **Separation**: the two exchanged only a sanitized status file (steps, pass/fail, error classes; coordinates,
  genotypes, HGVS, HPO terms, sample IDs and gene names masked).
- **Model use on patient data**: <none / list the Runner sessions and what they were for>.

## 8. Acknowledgement

<Copy the acknowledgement text required at the bottom of the Official Rules page.>
