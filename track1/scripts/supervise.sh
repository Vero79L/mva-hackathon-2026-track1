#!/usr/bin/env bash
# supervise.sh — Vero 在獨立的終端機視窗自己執行（Builder / Claude Code 不得執行本腳本，
# 也不得讀 work/ 裡的 log）。依序跑 03 → 事先登錄檢查 → 04 → 05 → 06 → 07 → 08；
# 每一步的完整輸出照常顯示在終端機；每步結束把遮罩後的狀態寫進 ~/mva-hackathon/handoff/STATUS.md。
# 失敗就分類（CLAUDE.md 第 5 節的 8 類）、寫入遮罩過的最後一行錯誤訊息，然後停下來。
# 重跑直接再執行一次即可：各步驟會自動跳過已完成的部分。
. "$(dirname "$0")/lib.sh"
. "$SCRIPTS/mask.sh"
stay_awake "$@"

HANDOFF="$MVA_ROOT/handoff"
STATUS="$HANDOFF/STATUS.md"
mkdir -p "$HANDOFF" "$WORK/logs"
[ -d "$RES" ] || die "參考資源資料夾 $RES 不存在：先接上外接碟（DISK-2）再跑。"

TABLE=""   # bash 3.2：用累積字串存表格列

write_status() {   # write_status <已遮罩的 message｜-> <next>
  {
    printf '# Track 1 status (sanitized)\n'
    printf 'updated: %s\n' "$(date '+%Y-%m-%d %H:%M')"
    printf '| step | result | started | minutes | error class | where |\n'
    printf '|---|---|---|---|---|---|\n'
    printf '%s' "$TABLE"
    if [ "$1" != "-" ]; then printf 'message (masked): %s\n' "$1"; fi
    printf 'next: %s\n' "$2"
  } > "$STATUS"
}

classify() {   # classify <log 檔> → 8 類之一
  if grep -qiE 'command not found|找不到指令|No module named|ModuleNotFoundError|不存在或無法啟動' "$1"; then echo tool-missing
  elif grep -qiE 'curl: \([0-9]+\)|下載失敗|下載不完整|Could not resolve|Connection (refused|reset|timed out)|transfer closed' "$1"; then echo download
  elif grep -qiE 'No space left|磁碟只剩|Disk quota' "$1"; then echo disk-space
  elif grep -qiE 'OutOfMemoryError|Cannot allocate memory|Killed: 9|heap space' "$1"; then echo memory
  elif grep -qiE 'Permission denied|Operation not permitted|Read-only file system' "$1"; then echo permission
  elif grep -qiE 'unsorted|not (sorted|indexed)|malformed|parse error|invalid (header|format|record)|名稱、長度|不一致' "$1"; then echo input-format
  elif grep -qiE 'syntax error|unbound variable|Traceback|AssertionError|KeyError|IndexError' "$1"; then echo script-bug
  else echo unknown; fi
}

run_step() {   # run_step <步驟名>（對應 scripts/<步驟名>.sh）
  local s="$1" start t0 mins rc log cls msg
  start="$(date '+%H:%M')"
  t0=$SECONDS
  log="$WORK/logs/supervise_${s}.last.out"
  say "supervisor：開始 ${s}"
  set +e
  bash "$SCRIPTS/$s.sh" 2>&1 | tee "$log"
  rc=${PIPESTATUS[0]}
  set -e
  mins=$(( (SECONDS - t0 + 59) / 60 ))
  if [ "$rc" = 0 ]; then
    TABLE="${TABLE}| $s | ok | $start | $mins | | |
"
    write_status "-" "running"
  else
    cls="$(classify "$log")"
    msg="$(tail -20 "$log" | grep -vE '^[[:space:]]*$' | tail -1 | mva_mask)"
    TABLE="${TABLE}| $s | failed | $start | $mins | $cls | $s.sh |
"
    write_status "$msg" "waiting for Builder fix"
    die "步驟 $s 失敗（class: ${cls}）。狀態已寫進 ${STATUS}；跟 Builder 說「看 STATUS.md」，修好後重跑本腳本。"
  fi
}

say "supervisor 開始（狀態檔：${STATUS}）"
run_step 03_fetch_data
mva_mask_init   # 03 之後才有 VCF：收集樣本名與 gated-data 檔名供遮罩

# 事先登錄關卡：規則沒 commit 就停（03 之後、04 之前）
if ! ( check_preregistered ); then
  TABLE="${TABLE}| preregistration | failed | $(date '+%H:%M') | 0 | script-bug | lib.sh |
"
  write_status "ranking_rules.yml not committed (or has uncommitted changes)" "waiting for Builder to commit track1"
  die "排序規則還沒 commit（或有未 commit 的修改）。讓 Builder commit track1 之後重跑本腳本。"
fi

for s in 04_qc 05_annotate 06_exomiser 07_rank 08_submission; do
  run_step "$s"
done

write_status "-" "all steps done — review the draft CSV and self-score in results/track1, then upload"
say "全部完成。提交草稿與自我評分在 ${RESULTS}；確認沒問題就由 Vero 上傳 CSV。"
