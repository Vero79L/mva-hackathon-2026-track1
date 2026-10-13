#!/usr/bin/env bash
# 10 Phase 2 analyses and the methods report. Reads patient data: Vero runs it in her own terminal
# (never from Claude Code). Order:
#   1. config/phase2_protocol.yml must be committed (pre-registration of the two analyses)
#   2. sensitivity.py      re-runs route B with one rule changed at a time (does not touch $WORK outputs)
#   3. aneuploidy_scan.py  per-chromosome depth and allele balance from norm.vcf.gz
#   4. report/fill_report.py  fills the report and exports the PDF (Chrome, Chromium, Edge or Brave)
# Output: results/track1/report/<team>_track1_report.html and .pdf
# Usage: bash track1/scripts/10_report.sh [--github URL] [--author "Name"] [--ai-training off|on]
#        [--lb-submission 1 --lb-rank-points 100 --lb-fmax 1.0 --lb-match full]
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
start_log 10_report
conda_on mva
need bcftools "（先跑 01_setup_mac.sh）"

say "1/4 事先登錄檢查（規則和分析計畫都要已 commit）"
check_preregistered
proto="$CONF/phase2_protocol.yml"
[ -f "$proto" ] || die "缺 config/phase2_protocol.yml"
pc="$(git -C "$CONF" log -1 --format='%h %ci' -- phase2_protocol.yml 2>/dev/null || true)"
[ -n "$pc" ] || die "phase2_protocol.yml 還沒 commit。先在終端機執行：
    cd \"$REPO\" && git add track1 && git commit -m \"Track 1: pre-register phase 2 analyses\"
  然後再跑一次。"
git -C "$CONF" diff --quiet HEAD -- phase2_protocol.yml || die "phase2_protocol.yml 有還沒 commit 的修改：先 commit（改計畫請把 version +1）"
ok "分析計畫已事先登錄：commit ${pc}"
[ -s "$WORK/trackB_pairs.tsv" ] || die "缺 trackB_pairs.tsv：先跑 supervise.sh（或 run_phase1.sh）"
[ -s "$WORK/norm.vcf.gz" ] || die "缺 norm.vcf.gz：先跑 05_annotate.sh"
out="$RESULTS/report"
mkdir -p "$out"

say "2/4 敏感度分析（約 1–3 分鐘）"
rc=0
python "$SCRIPTS/sensitivity.py" --work "$WORK" --rules "$CONF/ranking_rules.yml" --protocol "$proto" --out "$out/sensitivity.json" || rc=$?
if [ "$rc" -eq 3 ]; then
  warn "重跑凍結規則沒有重現原本的第 1 名：報告會標示敏感度結果不能解讀。先確認 work 裡的表在 07 之後沒被改過。"
elif [ "$rc" -ne 0 ]; then
  die "sensitivity.py 失敗（exit ${rc}），看上面的錯誤訊息"
fi
manifest 10_sensitivity "protocol ${pc%% *}"

say "3/4 非整倍體掃描（約 3–10 分鐘）"
sex="$(cat "$WORK/sex.txt" 2>/dev/null || echo UNKNOWN)"
python "$SCRIPTS/aneuploidy_scan.py" --vcf "$WORK/norm.vcf.gz" --sex "$sex" --protocol "$proto" --out "$out/aneuploidy.json"
manifest 10_aneuploidy "protocol ${pc%% *}"

[ -s "$T1/benchmark/output/benchmark.json" ] \
  || warn "還沒有基準測試結果：先跑 09_benchmark.sh（只用公開資料，約 1 小時），報告的「可擴展性」一節才會有數字"

say "4/4 填入報告、主辦方的 methods 表格、輸出 PDF"
python "$T1/report/fill_report.py" --work "$WORK" --results "$RESULTS" --t1 "$T1" --out-dir "$out" \
  --proband "$PROBAND" --hpo-obo "$HPO_OBO" "$@"
say "完成。報告在 ${out}：先打開 HTML 檢查，再上傳 PDF（或 HTML 存成 PDF）。
    填好的 methods 表格也寫進了 track1/docs/（不含基因型資料），檢查後跟著 track1 一起 commit。"
