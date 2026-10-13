# Track 1 — Variant prioritization for PROBAND01 (Phase 1)

Pipeline for **Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1**: rank candidate causal variants for
PROBAND01 from the organizers' WGS VCF (GRCh38) and produce a submission CSV of at most 10 rows.

It runs two independent routes and merges them only at the end:

| Route | What it does | Why |
|---|---|---|
| **A — Exomiser 15.1.0** (data 2602) | Genome-wide, phenotype-driven (HPO) ranking, run twice: the `exome` preset and the exome-with-introns preset | Makes no assumption about the gene |
| **B — pre-registered rules** (`config/ranking_rules.yml`) | Scores every variant in 13 MVA / mitotic-checkpoint genes (±10 kb, coding and non-coding, PASS and non-PASS) with a simplified ACMG/AMP point system, adds literature hypomorphic alleles that sit outside that window (e.g. the BUB1B upstream allele NC_000015.10:g.40117088G>A, about 44 kb upstream), then pairs alleles within each gene and excludes pairs proven to be in *cis* by read-backed phasing (PGT/PID) | Exomiser's exome presets only analyse PASS calls and drop upstream and deep-intronic variants, which is exactly where a hypomorphic second MVA allele can sit |

**Pre-registration.** The ranking rules are committed to git before any result is viewed. `07_rank.sh` and
`run_phase1.sh` refuse to run if `ranking_rules.yml` is uncommitted or has uncommitted edits. Every ranking run
is logged with the rules commit hash, version and SHA-256 (`results/track1/run_manifest.tsv`), so the report can
show that the rules predate the results.

## Data boundary

No patient data is stored in this repository. Raw data, intermediate files with per-site genotypes and logs live
outside the repo in `~/mva-hackathon/{gated-data,work}` and are deleted after the competition. Submission files and
phenotype lists live in `~/mva-hackathon/results`. The repo holds only code, rules, software versions (`env/`) and
genome-wide aggregate QC numbers (`qc/`). `.gitignore` and an optional pre-commit hook (`scripts/git-pre-commit.sh`)
block genomic files, genotype tables and files over 5 MB.

LLM use follows the organizers' guidance (Community #2). Code and documentation were written with Claude on a
subscription plan (claude.ai and Claude Code). Neither had access to patient data: the chat was connected only to the
code folders, and Claude Code ran with read-deny rules and an OS-level sandbox over the data folders. Steps that read
data run locally as plain scripts; when one fails, an agent built on the Claude Agent SDK and authenticated with an
API key (Anthropic Commercial Terms) diagnoses environment problems on the same machine. The two exchange only a
sanitized status file. See `CLAUDE.md`.

## Layout

```
~/mva-hackathon/
├── gated-data/            competition data (not in git; deleted after the competition)
├── work/track1/           intermediate files with genotypes (not in git; deleted after the competition)
├── resources/             public references: FASTA, VEP cache, ClinVar, Exomiser, hp.obo, official scorer
├── results/track1/        submission CSVs, HPO list, gene-level Exomiser summary, run manifest
└── repo/track1/           ← this folder
    ├── config/            paths.env, ranking_rules.yml (pre-registered), splitvep_format.txt
    ├── scripts/           00–08, run_phase1.sh, Python tools
    ├── tests/             synthetic-data tests (no competition data needed)
    ├── docs/              methods report draft
    ├── env/  qc/          versions and aggregate QC, written by the scripts
```

## Run order (macOS, bash 3.2 compatible; steps 02–08 also run on Linux)

| Step | Command | What | Time |
|---|---|---|---|
| 0 | `bash track1/scripts/00_check.sh` | Hardware, FileVault, Time Machine, tools, resources, progress. Read-only; run any time | 1 min |
| 1 | `bash track1/scripts/01_setup_mac.sh` | Xcode CLT, Homebrew, Java 21, Miniforge, conda envs `mva` / `vep` / `spliceai`, version records, Time Machine exclusions, git hook | 30–60 min |
| 2 | `bash track1/scripts/02_fetch_resources.sh` | Reference FASTA (chr prefix removed), ClinVar, panel BED, hp.obo, official scorer, Exomiser + data, VEP cache | overnight |
| 3 | `bash track1/scripts/03_fetch_data.sh` | Dataset manifest, selective download (no FASTQ/BAM), VCF header checks, HPO extraction | 10 min |
| — | `git add track1 && git commit` | **Pre-register the rules** before step 4 | — |
| 4–8 | `bash track1/scripts/run_phase1.sh` | QC + sex → normalization + VEP + SpliceAI + SV overlap → Exomiser ×2 → route B ranking → draft CSV + format check + self-scoring | 1.5–3 h |

Steps 4–8 can also be run one at a time (`04_qc.sh` … `08_submission.sh`). Finished steps whose inputs have not
changed are skipped; `FORCE=1` reruns everything. Long steps keep the Mac awake with `caffeinate`.

## Configuration

`config/paths.env` holds every path and version. Personal overrides (memory, an external SSD, the VCF path)
go in `config/local.env`, which `00_check.sh` creates from the machine's RAM and which is not committed. Anything
derived from the child's data stays outside the repo; a manual sex call, if ever needed, goes in
`~/mva-hackathon/work/track1/sex_override.txt`.

## Tests

`bash track1/tests/run_tests.sh` runs 31 tests on synthetic data in about 10 seconds: rule scoring and tie-breaks,
cis exclusion, multi-allelic splits, the known-hypomorph rule, SpliceAI gene matching, gene aliases, CSV assembly and
validation, scorer simulation, the pre-registration gate, downloads (resuming, a new version on the server, reruns
while offline), and macOS compatibility of every shell script (including `${var}` before Chinese text, which the
built-in bash 3.2 misreads).

---

## 中文快速說明

1. 先跑 `00_check.sh` 看哪裡還沒好；它會依記憶體自動建立 `config/local.env`。
2. `01_setup_mac.sh` 裝工具（會要你輸入 Mac 密碼）；`02_fetch_resources.sh` 晚上接電源掛著下載。
3. `03_fetch_data.sh`：已經有 VCF 就不會重抓；HPO 清單放在 `results/track1/hpo_terms.tsv`，不進 repo。
4. **跑 04 之前先 commit 規則**：`cd ~/mva-hackathon/repo && git add track1 && git commit -m "Track 1: pipeline + pre-registered ranking rules v1"`。
5. `run_phase1.sh` 一次跑完 04–08。最後畫面會列出 10 列草稿、格式檢查、和「如果答案是第 k 列會拿幾分」的模擬。
6. 提交前人工看過第 1 列的兩個等位基因（`work/track1/trackB_variants.tsv`）和 `panel_sv.tsv`。
7. 要改規則：`version` +1、另外 commit，報告寫明改了什麼、為什麼改。
8. 讓 Claude 在 Mac 上幫忙時分兩個角色（細節在 `CLAUDE.md`）：
   - **Builder**：Claude Code，用 Max 訂閱登入。只寫程式、跑 00–02 和測試，只讀 `~/mva-hackathon/handoff/STATUS.md`，碰不到資料（用使用者層級的 deny 規則和 sandbox 擋）。
   - **Runner**：妳在另一個終端機視窗啟動 supervisor，由它直接跑 03 → commit → `run_phase1.sh`；某一步失敗才叫 Claude Agent SDK 進來診斷。SDK 用 Keychain 裡的 API key（API 條款，扣訂閱附送的 API credits），只修環境、不改程式和規則。

完整說明見 Phase 1 執行 Spec v1.2（filing 資料夾裡的 `specs/phase-1/phase1-spec-v1.2.html`）。
