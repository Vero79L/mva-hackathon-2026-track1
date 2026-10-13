# Track 1 runbook: Builder (Claude Code) and Runner (Agent SDK) on Vero's Mac

You are helping Vero run Phase 1 of the Rare Disease, Real Kid: MVA Hackathon 2026, Track 1, on her own Mac.
Reply to Vero in Traditional Chinese. The steps themselves are defined in Phase 1 Spec v1.2 (`README.md` here, and
the spec HTML in her filing folder). This file says who runs which step and where the data boundary is.
Where this file and the spec differ on who runs what, this file wins.

## 0. Which role are you?

Two Claude processes work in this folder. Work out which one you are before doing anything.

| | Builder | Runner |
|---|---|---|
| What it is | The interactive Claude Code session | A Claude Agent SDK agent started by the supervisor |
| Account | Vero's Max subscription (consumer terms) | API key from Vero's Console workspace (API terms, paid from the subscriber API credits) |
| How to tell | `/status` shows a claude.ai login | The supervisor sets `MVA_ROLE=runner` and says so in your system prompt |
| Reads patient data | Never | Yes, on this Mac only |
| Writes code | Yes: `track1/scripts`, `track1/tests`, the supervisor | No |
| Runs | Steps 0, 1, 1a, 2 and the tests | Diagnoses failed steps 3–8 |

If you can't tell, act as the Builder.

The competition allows patient data to reach an LLM only under business / API terms (organizers' answer in
Community #2; #21 on consumer plans with training off is still open). That is why only the Runner may see data.

## 1. Data boundary

Patient data is everything under `~/mva-hackathon/gated-data`, `~/mva-hackathon/work` and `~/mva-hackathon/results`:
the VCF, logs, Exomiser output, the HPO list, the phenopacket, the submission CSV. Treat `track1/qc/` the same way:
its numbers come from the child's genome.

- **Builder**: never open, grep, list contents of, or run a data step (03–08) on those folders. The only thing you read
  from the pipeline is `~/mva-hackathon/handoff/STATUS.md`. Vero's user-level settings enforce this (section 6); don't
  ask her to loosen them.
- **Runner**: don't use web tools, don't upload or paste data anywhere, and don't write patient-derived values into
  the repo (including `track1/config/local.env`). If the sex check is UNKNOWN and Vero gives you the answer, write it to
  `~/mva-hackathon/work/track1/sex_override.txt` (MALE or FEMALE).
- **Both**: the only channel from Runner to Builder is `STATUS.md` in the format of section 5.

## 2. Who runs what

Every script is safe to rerun: downloads resume, finished steps whose inputs did not change are skipped, `FORCE=1`
redoes everything. Each script logs to `~/mva-hackathon/work/track1/logs/`. Paths are full, so they work from anywhere.

| Step | Command | Reads data | Who | Typical time |
|---|---|---|---|---|
| 0 | `bash ~/mva-hackathon/repo/track1/scripts/00_check.sh` | no | Builder | 1 min |
| 1 | `bash ~/mva-hackathon/repo/track1/scripts/01_setup_mac.sh` | no | Builder (Vero types passwords) | 30–60 min |
| 1a | `bash ~/mva-hackathon/repo/track1/tests/run_tests.sh` (expect `OK`; 35 tests) | no | Builder | 10 s |
| 2 | `bash ~/mva-hackathon/repo/track1/scripts/02_fetch_resources.sh` | no | Builder | hours |
| 3 | `bash ~/mva-hackathon/repo/track1/scripts/03_fetch_data.sh` | yes | supervisor | 10 min |
| 3a | pre-register the rules (section 3) | no | Builder or Vero | 1 min |
| 4–8 | `bash ~/mva-hackathon/repo/track1/scripts/run_phase1.sh` | yes | supervisor | 1.5–3 h |
| 9 | upload CSV + report | — | Vero | 15 min |

The supervisor:

- Vero starts it herself in a separate terminal window. The Builder must not start it: anything the Builder starts
  sends its output back into the Builder's conversation, and inside the Builder's sandbox it can't read the data.
- It runs the scripts directly, without an LLM. It starts the Runner only when a step exits non-zero, and pauses
  after step 3 until track1 is committed (3a).
- It reads the API key from the macOS Keychain (`security find-generic-password -a vero -s mva-anthropic-api -w`) and
  passes it only to the Runner process. Never write the key to a file, a log, or `~/.zshrc`.
- Runner process environment: `MVA_ROLE=runner`,
  `CLAUDE_CONFIG_DIR=~/mva-hackathon/work/track1/.agent-config` (keeps the Runner's transcripts and memory out of
  the Builder's `~/.claude`), `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1`. Options: `setting_sources=["project"]`, Sonnet,
  no WebFetch / WebSearch, a `max_turns` cap, deny rules from section 3.

## 3. Pre-registration (do not skip, do not bend)

- `track1/config/ranking_rules.yml` must be committed before step 4; `07_rank.sh` and `run_phase1.sh` refuse to run
  otherwise. Commit only the track1 folder:
  `git -C ~/mva-hackathon/repo add track1 && git -C ~/mva-hackathon/repo commit -m "Track 1: pipeline + pre-registered ranking rules v1"`
- The Runner sees data, so it must not change ranking or submission logic. Deny it Edit/Write on
  `track1/scripts/**`, `track1/tests/**`, `track1/config/ranking_rules.yml`, and deny `git commit` / `git push`.
- Code bugs found during steps 3–8 go into STATUS.md. The Builder, who has not seen the data, fixes the code, runs the
  tests and commits; then the supervisor reruns. This keeps the claim "code and rules were written before looking at
  the data" true for the methods report.
- After the first ranking run, changing `ranking_rules.yml` needs Vero's explicit approval, `version: 2`, a separate
  commit, and a written reason for the report.

## 4. Git safety and failures

- Never run `git add .` or `git add -A` at the repo root: `config/cases/` and `evidence/` belong to Track 2 and hold
  case-derived information.
- Never copy anything from `gated-data/`, `work/` or `results/` into the repo; don't bypass the pre-commit hook.
- Don't push to a public repository and don't submit to the competition. Vero uploads the CSV.
- Runner fixes only the environment: missing tools, PATH, downloads, disk space (`RES` in `local.env`), memory
  (`JAVA_XMX`, `VEP_FORK` in `local.env`). Anything else → STATUS.md.
- Builder fixes code with the smallest change, keeps macOS bash 3.2 / BSD compatibility (no `sed -i`, `grep -P`,
  `readlink -f`, `declare -A`, `| grep -q` under pipefail; write `${var}` whenever Chinese text follows a variable,
  because bash 3.2 reads the first bytes of the character as part of the name), reruns
  `bash ~/mva-hackathon/repo/track1/tests/run_tests.sh`, commits with a `fix:` message.
- Known fallbacks: no Apple Silicon build in bioconda → `01_setup_mac.sh` retries with `CONDA_SUBDIR=osx-64` (needs
  Rosetta); no `ensembl-vep=116` → set `VEP_VERSION` in `local.env` to the newest available, rerun 01 and 02;
  Exomiser `OutOfMemoryError` → lower `JAVA_XMX`; SpliceAI won't install → continue, it is optional.

## 5. STATUS.md: the only handoff

Path: `~/mva-hackathon/handoff/STATUS.md` (outside the three data folders). The supervisor's code writes it from
fixed fields, then masks the free-text message; don't rely on an LLM alone to sanitize.

```text
# Track 1 status (sanitized)
updated: <local time>
| step | result | started | minutes | error class | where |
|---|---|---|---|---|---|
| 03_fetch_data | ok | 14:20 | 6 | | |
| 04_qc | failed | 14:27 | 2 | tool-missing | 04_qc.sh:31 |
message (masked): bcftools: command not found
next: waiting for Builder fix
```

- Error classes: `tool-missing`, `download`, `disk-space`, `memory`, `permission`, `input-format`, `script-bug`,
  `unknown`.
- Mask before writing: genomic coordinates (`15:123`, `chr15:123`), any run of 4 or more digits, alleles
  (`A>G`, runs of ACGTN), HGVS (`c.`, `p.`, `g.`, `n.`), rsIDs, `HP:` IDs, genotypes (`0/1`, `1|0`), sample names
  (read them from the VCF header at run time), the 13 MVA gene names, `MALE` / `FEMALE`, and file names under
  `gated-data`.
- Vero may paste STATUS.md into her claude.ai project chat. Leaderboard scores after submission may be added.

## 6. Builder settings (Vero's user-level `~/.claude/settings.json`)

Put these in the user-level file, not in a repo `.claude/settings.json`: the Runner loads project settings and would
be locked out of the data too.

```json
{
  "permissions": {
    "deny": [
      "Read(~/mva-hackathon/gated-data/**)", "Edit(~/mva-hackathon/gated-data/**)",
      "Read(~/mva-hackathon/work/**)", "Edit(~/mva-hackathon/work/**)",
      "Read(~/mva-hackathon/results/**)", "Edit(~/mva-hackathon/results/**)",
      "Read(~/mva-hackathon/repo/track1/qc/**)"
    ]
  },
  "sandbox": {
    "enabled": true,
    "filesystem": {
      "denyRead": ["~/mva-hackathon/gated-data", "~/mva-hackathon/work", "~/mva-hackathon/results"]
    }
  }
}
```

The deny rules cover Claude Code's file tools; the sandbox covers shell commands. If the sandbox asks to run a
command outside it, Vero approves only commands that don't touch those folders.
