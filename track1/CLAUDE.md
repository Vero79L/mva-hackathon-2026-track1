# Track 1 runbook for Claude Code (local, on Vero's Mac)

You are helping Vero run Phase 1 of the Rare Disease, Real Kid: MVA Hackathon 2026, Track 1, on her own Mac.
Reply to Vero in Traditional Chinese. The full plan is Phase 1 Spec v1.2 (`README.md` here, and the spec HTML in her
filing folder); this file is the short operating version.

## 0. Check the account before touching any data

The competition's data rules allow patient data to be processed by an LLM only under business / API terms
(the organizers' answer in Community discussion #2). A Claude Pro or Max login does not qualify.

1. Run `/status`. It must show a Claude Console (API) credential: an Anthropic profile or an API key.
2. If it shows a claude.ai subscription login, you may run only the steps that do not read the child's data
   (steps 0, 1, 1a and 2 in the table below). Stop before step 03 and ask Vero
   to switch with `/login` → Anthropic Console account.
3. Never send the child's data anywhere else (no web tools, no uploads, no pastes into other services).

## 1. Run order

Vero starts Claude Code in `~/mva-hackathon/repo/track1` so this file loads; the commands below use full paths, so
they work from any directory. Every script is safe to rerun: downloads resume, finished steps whose inputs did not
change are skipped, `FORCE=1` redoes everything. Each script logs to `~/mva-hackathon/work/track1/logs/`.

| Step | Command | Reads patient data | Typical time |
|---|---|---|---|
| 0 | `bash ~/mva-hackathon/repo/track1/scripts/00_check.sh` | no | 1 min |
| 1 | `bash ~/mva-hackathon/repo/track1/scripts/01_setup_mac.sh` | no | 30–60 min |
| 1a | `bash ~/mva-hackathon/repo/track1/tests/run_tests.sh` (expect 30 tests OK) | no | 10 s |
| 2 | `bash ~/mva-hackathon/repo/track1/scripts/02_fetch_resources.sh` | no | hours (VEP cache, Exomiser data) |
| 3 | `bash ~/mva-hackathon/repo/track1/scripts/03_fetch_data.sh` | yes | 10 min |
| 3a | pre-register the rules (see section 2) | no | 1 min |
| 4–8 | `bash ~/mva-hackathon/repo/track1/scripts/run_phase1.sh` | yes | 1.5–3 h |

- Long steps: run them in the background and check the newest log in `work/track1/logs/` instead of blocking.
- `01_setup_mac.sh` asks for Vero's Mac password (Homebrew, Java). Ask her to type it in her own terminal. Never
  ask for the password in chat.
- `hf auth login` (Hugging Face token) is also something Vero types herself.
- `02_fetch_resources.sh` needs about 150 GB free. If space is short, set `RES` in `track1/config/local.env`
  (see `config/paths.env`), don't delete anything.

## 2. Pre-registration (do not skip, do not bend)

- `track1/config/ranking_rules.yml` must be committed before step 4. `07_rank.sh` and `run_phase1.sh` refuse to run
  otherwise. Commit only the track1 folder:
  `git -C ~/mva-hackathon/repo add track1 && git -C ~/mva-hackathon/repo commit -m "Track 1: pipeline + pre-registered ranking rules v1"`
- After the first ranking run, do not edit `ranking_rules.yml` to change which variants come out on top. A rules
  change needs Vero's explicit approval, `version: 2`, a separate commit, and a written reason for the report.

## 3. Git and data safety

- Never run `git add .` or `git add -A` at the repo root: `config/cases/` and `evidence/` belong to Track 2 and
  hold case-derived information.
- Never copy anything from `gated-data/`, `work/` or `results/` into the repo. The pre-commit hook and
  `track1/.gitignore` block genomic files, genotype tables and files over 5 MB; don't bypass them with `--no-verify`
  for data.
- Don't push to a public repository and don't submit to the competition. Prepare the CSV; Vero uploads it.

## 4. When something fails

- Read the error and the newest log, fix the cause, rerun the same script.
- Script bugs: make the smallest fix in `track1/scripts`, keep macOS bash 3.2 / BSD compatibility (no `sed -i`,
  `grep -P`, `readlink -f`, `declare -A`, `| grep -q` under pipefail), then rerun `bash ~/mva-hackathon/repo/track1/tests/run_tests.sh`.
- Known environment fallbacks: no Apple Silicon build in bioconda → `01_setup_mac.sh` retries with `CONDA_SUBDIR=osx-64`
  (needs Rosetta); no `ensembl-vep=116` → set `VEP_VERSION` in `local.env` to the newest available and rerun 01 and 02;
  Exomiser `OutOfMemoryError` → lower `JAVA_XMX` in `local.env`; SpliceAI fails to install → continue, it is optional.

## 5. Report back

After each step, tell Vero in Traditional Chinese what ran, how long it took, and what is next. Vero also keeps a
project log in a claude.ai conversation, which must not receive patient data. When she asks for a status to paste
there, write a block titled `Status for claude.ai` that contains only: steps run, pass / fail, tool versions,
run times, and error messages with every patient-derived value removed. Leave out variant coordinates, genotypes,
HGVS, gene-level results for this child, HPO terms, sample IDs, QC numbers and leaderboard-relevant details of the
candidates. Scores from the leaderboard after submission may be included.
