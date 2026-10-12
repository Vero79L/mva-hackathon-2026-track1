#!/usr/bin/env bash
# 第一階段「一行指令重跑」：04 QC → 05 標準化 / 註解 → 06 Exomiser → 07 路線 B 排序 → 08 提交草稿與自我評分
# 前提：01–03 已完成（工具、參考資源、比賽資料），ranking_rules.yml 已 commit。
# 已完成而且輸入沒變的步驟會自動跳過；FORCE=1 全部重做。全新跑一次約 1.5–3 小時（Exomiser 佔大半）。
# 用法：bash track1/scripts/run_phase1.sh
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
check_preregistered     # 先確認規則已事先登錄，避免跑了一小時才被 07 擋下
start=$(date +%s)
for s in 04_qc 05_annotate 06_exomiser 07_rank 08_submission; do
  bash "$SCRIPTS/$s.sh" || die "$s 失敗：修好後重跑 run_phase1.sh（已完成的步驟會自動跳過）"
done
say "第一階段全部完成，用了 $(( ($(date +%s) - start) / 60 )) 分鐘。提交前請人工看過第 1 列的兩個等位基因。"
