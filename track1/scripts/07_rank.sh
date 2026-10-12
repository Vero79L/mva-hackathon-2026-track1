#!/usr/bin/env bash
# 07 路線 B：依事先登錄的 ranking_rules.yml，對 MVA 基因組打分、配對、排序。
# 規則沒 commit（或 commit 後又改過）就不給跑：commit 時間早於第一次排序，是「規則先於結果」的證據。
# 每次執行都記進 $RESULTS/run_manifest.tsv（時間、規則 commit、版本、雜湊），寫方法報告時引用。
# 用法：bash track1/scripts/07_rank.sh
. "$(dirname "$0")/lib.sh"
start_log 07_rank
conda_on mva
check_preregistered
[ -s "$WORK/panel.vep.tsv" ] || die "缺 panel.vep.tsv：先跑 05_annotate.sh"

first="$RESULTS/first_rank_run.txt"
sha="$(sha256 "$CONF/ranking_rules.yml")"
if [ -s "$first" ] && ! grep "$sha" "$first" >/dev/null; then
  warn "規則和第一次排序時不同（$(head -1 "$first")）。這是規則 v$(rules_version)：方法報告要寫明改了什麼、為什麼改。"
fi

python "$SCRIPTS/rank_panel.py" --work "$WORK" --rules "$CONF/ranking_rules.yml"
manifest 07_rank "rank_panel.py"
if [ ! -s "$first" ]; then
  mkdir -p "$RESULTS"
  echo "first ranking run $(date '+%Y-%m-%d %H:%M:%S'); rules commit $RULES_COMMIT; rules v$(rules_version) sha256 $sha" > "$first"
  ok "第一次排序已記錄：$(cat "$first")"
fi
info "人工檢查前 5 名配對的每個等位基因：品質（FILTER、GQ、DP、AD）、HGVS、ClinVar、SpliceAI → $WORK/trackB_variants.tsv"
[ -s "$WORK/panel_sv.tsv" ] && warn "別忘了 SV / CNV：$WORK/panel_sv.tsv"
say "完成。下一步：bash track1/scripts/08_submission.sh"
