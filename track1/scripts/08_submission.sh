#!/usr/bin/env bash
# 08 組提交 CSV（≤10 列、EPCR 階梯、chr 前綴）→ 格式檢查 → 座標回查 VCF → 本機模擬與官方 evaluation.py 交叉檢查。
# 檔名：$RESULTS/<HF 帳號>_mva-tiered-v<規則版本>.csv（HF 帳號自動讀 hf 登入的帳號，或 HF_USER=xxx 指定）
# secondary findings：人工審核後填進 $RESULTS/secondary_findings.tsv（欄位：chrom pos ref alt gene reason），會放在第 9–10 列。
# 用法：bash track1/scripts/08_submission.sh            草稿 + 檢查
#       bash track1/scripts/08_submission.sh check <CSV>  只檢查一個現成的 CSV（例如手動改過之後）
. "$(dirname "$0")/lib.sh"
start_log 08_submission
conda_on mva
[ -s "$SCORER" ] || warn "沒有官方 evaluation.py（02_fetch_resources.sh scorer）：只做內建檢查"

if [ "${1:-}" = check ]; then
  csv="${2:?用法：08_submission.sh check <CSV 路徑>}"
else
  [ -s "$WORK/trackB_pairs.tsv" ] || die "缺 trackB_pairs.tsv：先跑 07_rank.sh"
  say "1/2 組提交草稿（依 ranking_rules.yml 的 merge 與 EPCR）"
  out="$(python "$SCRIPTS/make_submission.py" draft --work "$WORK" --results "$RESULTS" --rules "$CONF/ranking_rules.yml" \
         --proband "$PROBAND" | tee /dev/stderr | sed -nE 's/^CSV=(.*)$/\1/p')"
  csv="$out"
  [ -s "$csv" ] || die "沒有產生 CSV"
  manifest 08_draft "$(basename "$csv")"
fi

say "2/2 檢查 $csv"
python "$SCRIPTS/make_submission.py" check --csv "$csv" --vcf "$WORK/norm.vcf.gz" --scorer "$SCORER" --proband "$PROBAND"
say "提交前：每一列 notes 讀一遍、第 1 列的兩個等位基因回到 trackB_variants.tsv 再看一次品質。
    上傳：比賽 Space →「Submit - Track 1」→ 填 GitHub URL → 上傳 CSV 和方法報告 →「Submit & Score」。
    交完把結果（Rank points、F-max、match 類型）記進 Playbook 的提交紀錄。"
