# Track 1 — Variant prioritization for PROBAND01

Pipeline for **Rare Disease, Real Kid: MVA Hackathon 2026 — Track 1**: rank candidate causal variants for
PROBAND01 from the organizers' WGS VCF (GRCh38) and produce a submission CSV of at most 10 rows.

It runs two independent routes and merges them only at the end:

| Route | What it does | Why |
|---|---|---|
| **A — Exomiser 15.1.0** (data 2602) | Genome-wide, phenotype-driven (HPO) ranking, run twice: the `exome` preset and the exome-with-introns preset | Makes no assumption about the gene |
| **B — pre-registered rules** (`config/ranking_rules.yml`) | Scores every variant in 13 MVA / mitotic-checkpoint genes (±10 kb, coding and non-coding, PASS and non-PASS) with a simplified ACMG/AMP point system, adds literature hypomorphic alleles that sit outside that window (e.g. the BUB1B upstream allele NC_000015.10:g.40117088G>A, about 44 kb upstream), then pairs alleles within each gene and excludes pairs proven to be in *cis* by read-backed phasing (PGT/PID) | Exomiser's exome presets only analyse PASS calls and drop upstream and deep-intronic variants, which is exactly where a hypomorphic second MVA allele can sit |

**Pre-registration.** The ranking rules are committed to git before any result is viewed. `07_rank.sh` and
`run_phase1.sh` refuse to run if `ranking_rules.yml` is uncommitted or has uncommitted edits. Every ranking run
is logged with the rules commit hash, version and SHA-256 (`docs/preregistration/run_manifest.tsv`,
`docs/preregistration/first_rank_run.txt`), so the report can show that the rules predate the results. The Phase 2
analyses and the generalisation benchmark are pre-registered the same way in `config/phase2_protocol.yml`: steps 09
and 10 refuse to run until it is committed.

**Other diseases.** The ranking engine (`scripts/rank_panel.py` and the weights in `config/ranking_rules.yml`) holds no
disease knowledge. The disease enters only through the gene tiers and the list of known alleles, and
`scripts/make_rules.py` swaps just those, for example from a Genomics England PanelApp panel:

```
python track1/scripts/make_rules.py --genelists track1/benchmark/genelists_panelapp.tsv --panel fanconi_bloom \
  --out /tmp/fanconi_rules.yml
```

Step 09 tests this on public genomes (see below).

## Data boundary

No patient data is stored in this repository. Raw data, intermediate files with per-site genotypes and logs live
outside the repo in `~/mva-hackathon/{gated-data,work}` and are deleted after the competition. Submission files and
phenotype lists live in `~/mva-hackathon/results`. The repo holds only code, rules, software versions (`env/`) and
the pre-registration records (`docs/preregistration/`). The aggregate QC numbers in `qc/` come from the child's genome
and are git-ignored, like the filled report. `.gitignore` and an optional pre-commit hook (`scripts/git-pre-commit.sh`)
block genomic files, genotype tables and files over 5 MB.

LLM use follows the organizers' guidance (Community #2). Code and documentation were written with Claude on a
subscription plan (claude.ai and Claude Code). Claude Code ran with read-deny rules and an OS-level sandbox over the
data folders and never read patient data. Every step that reads data (03–08, then 10) is started by the participant
in a separate terminal (03–08 via `scripts/supervise.sh`), as plain bash scripts with no LLM component; on failure
supervise.sh writes a code-masked status file (`handoff/STATUS.md`) — the only channel back to the assistant — and
stops. Step 09 uses public data only. Full declaration:
`docs/Vero7878_track1_methods.md` §7. See also `CLAUDE.md`.

## Layout

```
~/mva-hackathon/
├── gated-data/            competition data (not in git; deleted after the competition)
├── work/track1/           intermediate files with genotypes (not in git; deleted after the competition)
├── resources/             public references: FASTA, VEP cache, ClinVar, Exomiser, hp.obo, official scorer
├── public-benchmark/      step 09 working files: GIAB genomes and their annotation (public data only)
├── results/track1/        submission CSVs, HPO list, gene-level Exomiser summary, run manifest, the filled report
└── repo/track1/           ← this folder
    ├── config/            paths.env, ranking_rules.yml and phase2_protocol.yml (pre-registered), splitvep_format.txt
    ├── scripts/           00–10, run_phase1.sh, supervise.sh, make_rules.py, Python tools
    ├── benchmark/         PanelApp gene lists (pinned versions) and the step 09 results in output/ (public data)
    ├── report/            report template and the script that fills it (step 10)
    ├── tests/             synthetic-data tests (no competition data needed)
    ├── docs/              methods file, the filled methods form, pre-registration records
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
| — | `git add track1 && git commit` | **Pre-register `phase2_protocol.yml`** (Phase 2 analyses and benchmark) before steps 9 and 10 | — |
| 9 | `bash track1/scripts/09_benchmark.sh` | Generalisation benchmark on public genomes only: GIAB HG001–HG007 (about 0.9 GB), ClinVar pathogenic genotypes placed into them, route B re-run under the MVA rules and seven PanelApp panels. Results in `benchmark/output/`; commit them | 1–2 h |
| 10 | `bash track1/scripts/10_report.sh --github <URL> --author "<name>" --ai-training <on\|off>` | Pre-registered sensitivity analysis (48 one-rule-at-a-time scenarios) and mosaic-aneuploidy scan, then the report (HTML with a navigator and interactive figures, plus a PDF printed by Chrome) and the organisers' methods form | 5–15 min |

Steps 4–8 can also be run one at a time (`04_qc.sh` … `08_submission.sh`). Finished steps whose inputs have not
changed are skipped; `FORCE=1` reruns everything. Long steps keep the Mac awake with `caffeinate`.

## Configuration

`config/paths.env` holds every path and version. Personal overrides (memory, an external SSD, the VCF path)
go in `config/local.env`, which `00_check.sh` creates from the machine's RAM and which is not committed. Anything
derived from the child's data stays outside the repo; a manual sex call, if ever needed, goes in
`~/mva-hackathon/work/track1/sex_override.txt`.

## Generalisation benchmark (step 09)

Does route B work for other people and other recessive diseases, with nothing changed but the gene list? Step 09
answers this on public data, as pre-registered in `config/phase2_protocol.yml` (`benchmark:`). It places ClinVar
pathogenic or likely pathogenic variants, as compound-heterozygous or homozygous genotypes, into seven GIAB genomes and
re-runs the frozen `rank_panel.py` under the MVA rules and under rules built by `make_rules.py` from seven PanelApp
panels (`benchmark/genelists_panelapp.tsv`, pinned versions, 1 to about 1,500 green genes). Each genotype is ranked
twice: with its ClinVar evidence (positive control) and with ClinVar removed, as for a newly found variant (the
primary test). The same genomes with nothing placed count the background candidates per panel. Annotation uses the
same VEP release and flags as step 05. Results go to `benchmark/output/benchmark.json` and `simulations.tsv`, public
data only, and are committed so the report can cite them.

## Report and methods form (step 10)

`report/report_template.html` holds the report text, the literature evidence and the figure code. `report/fill_report.py`
reads the run's outputs on the participant's machine, writes them into the template as one JSON block, and prints the
PDF with headless Chrome. Case-specific values therefore never pass through an assistant. The output,
`results/track1/report/<team>_track1_report.{html,pdf}`, is patient-derived and stays out of git. Besides the ranking,
the report covers evidence per allele, exclusions, the sensitivity analysis, the aneuploidy scan, scalability (the
benchmark, run time and cost), clinical translation and AI governance. The same script fills the organisers' methods
form (all 13 questions, abstract within 500 words) as `docs/<team>_track1_methods_form.{md,xlsx}`; it describes the
method and contains no genotype-level data.

## Tests

`bash track1/tests/run_tests.sh` runs 45 tests on synthetic data in about a minute: rule scoring and tie-breaks,
cis exclusion, multi-allelic splits, the known-hypomorph rule, SpliceAI gene matching, gene aliases, CSV assembly and
validation, scorer simulation, the pre-registration gate, downloads (resuming, a new version on the server, reruns
while offline), macOS compatibility of every shell script (including `${var}` before Chinese text, which the
built-in bash 3.2 misreads), the sensitivity analysis, the aneuploidy scan (40% and 15% mosaic gains are
found, a normal genome stays silent; skipped without bcftools), the report and methods-form filler (no phenotype terms
or sex-chromosome metrics in the output, only submitted variants listed, abstract within 500 words) and the benchmark
(swapping the gene list leaves the engine identical, spike eligibility, and the expected rank of each placed genotype
on synthetic genomes).

---

## 中文快速說明

1. 先跑 `00_check.sh` 看哪裡還沒好；它會依記憶體自動建立 `config/local.env`。
2. `01_setup_mac.sh` 裝工具（會要你輸入 Mac 密碼）；`02_fetch_resources.sh` 晚上接電源掛著下載。
3. `03_fetch_data.sh`：已經有 VCF 就不會重抓；HPO 清單放在 `results/track1/hpo_terms.tsv`，不進 repo。
4. **跑 04 之前先 commit 規則**：`cd ~/mva-hackathon/repo && git add track1 && git commit -m "Track 1: pipeline + pre-registered ranking rules v1"`。
5. `run_phase1.sh` 一次跑完 04–08。最後畫面會列出 10 列草稿、格式檢查、和「如果答案是第 k 列會拿幾分」的模擬。
6. 提交前人工看過第 1 列的兩個等位基因（`work/track1/trackB_variants.tsv`）和 `panel_sv.tsv`。
7. 要改規則：`version` +1、另外 commit，報告寫明改了什麼、為什麼改。
7a. 第二階段：先 commit `config/phase2_protocol.yml`（事先登錄），再跑 `09_benchmark.sh`（只用公開資料，約 1–2 小時，跑完 commit `track1/benchmark`），最後自己在終端機跑 `10_report.sh`（會讀資料，不要叫 Claude Code 跑）。報告在 `results/track1/report/`，不進 repo；填好的主辦方 methods 表格寫進 `docs/`，跟 track1 一起 commit。
8. 讓 Claude 在 Mac 上幫忙時只有一個角色（細節在 `CLAUDE.md`）：
   - **Builder**：Claude Code，用 Max 訂閱登入。只寫程式、跑 00–02 和測試，只讀 `~/mva-hackathon/handoff/STATUS.md`，碰不到資料（用使用者層級的 deny 規則和 sandbox 擋）。
   - 03–08 由妳自己在另一個終端機跑 `scripts/supervise.sh`（純腳本，沒有 LLM）；某一步失敗時它把遮罩後的狀態寫進 STATUS.md 就停，Builder 修好程式後由妳重跑。

完整說明見 Phase 1 執行 Spec v1.2（filing 資料夾裡的 `specs/phase-1/phase1-spec-v1.2.html`）。
