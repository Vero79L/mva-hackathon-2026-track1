#!/usr/bin/env bash
# 03 比賽資料：列出資料集清單 → 只在缺主 VCF 時選擇性下載（不抓 FASTQ / BAM / CRAM）→ 檢查 VCF 檔頭 → 整理 HPO terms。
# 已經下載過的資料（$DATA 裡已有主 VCF）不會重抓。
# 用法：bash track1/scripts/03_fetch_data.sh
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
start_log 03_fetch_data
conda_on mva
mkdir -p "$DATA" "$WORK" "$RESULTS"

say "1/4 Hugging Face 登入"
need hf "（mva 環境裡應該有 hf 指令：重跑 01_setup_mac.sh）"
hf auth whoami >/dev/null 2>&1 || die "還沒登入 Hugging Face：conda activate mva && hf auth login（貼 Read token；問要不要存成 git credential 選 n）"
ok "已登入"

say "2/4 資料集清單（只有檔名與大小）→ $RESULTS/dataset_manifest.txt"
python - "$HF_DATASET" > "$RESULTS/dataset_manifest.txt" <<'EOF'
import sys
from huggingface_hub import HfApi
info = HfApi().dataset_info(sys.argv[1], files_metadata=True)
for s in sorted(info.siblings, key=lambda x: x.rfilename):
    print(f"{(s.size or 0) / 1e9:8.2f} GB  {s.rfilename}")
EOF
cat "$RESULTS/dataset_manifest.txt"

say "3/4 主 VCF"
if vcf="$(find_vcf 2>/dev/null)"; then
  ok "$DATA 裡已經有主 VCF，不重新下載：${vcf#"$DATA"/}"
else
  info "選擇性下載：VCF、索引、表現型與說明文件（不含 FASTQ / BAM / CRAM）"
  need_space "$DATA" 5
  hf download "$HF_DATASET" --repo-type dataset --local-dir "$DATA" \
    --include "*.vcf.gz" "*.vcf.gz.tbi" "*.csi" "*.txt" "*.json" "*.tsv" "*.csv" "*.md" "*.pdf" "*.xlsx" "*.docx" "*.yml" "*.yaml" \
    --exclude "*.fastq*" "*.fq*" "*.bam" "*.bai" "*.cram" "*.crai"
  vcf="$(find_vcf)"
  ok "主 VCF：${vcf#"$DATA"/}"
fi
ensure_index "$vcf"
sv="$(sv_vcfs)"
if [ -n "$sv" ]; then info "SV / CNV VCF（05 會在 MVA 基因上比對）："; printf '%s\n' "$sv" | sed "s|^$DATA/|      |"; fi

say "4/4 VCF 檔頭檢查 → $WORK/vcf_header_check.txt"
hdr="$WORK/vcf_header_check.txt"
{
  echo "VCF: $vcf"
  bcftools view -h "$vcf" | grep -E '^##(fileformat|reference|source)' || true
  echo "samples: $(bcftools query -l "$vcf" | tr '\n' ' ')"
  echo "FORMAT: $(bcftools view -h "$vcf" | sed -nE 's/^##FORMAT=<ID=([^,]+).*/\1/p' | tr '\n' ' ')"
  echo "FILTER: $(bcftools view -h "$vcf" | sed -nE 's/^##FILTER=<ID=([^,]+).*/\1/p' | tr '\n' ' ')"
  echo "contigs (前 25): $(bcftools view -h "$vcf" | sed -nE 's/^##contig=<ID=([^,>]+).*/\1/p' | head -25 | tr '\n' ' ')"
} > "$hdr"
cat "$hdr"
ns="$(bcftools query -l "$vcf" | grep -c . || true)"
if [ "$ns" = 1 ]; then ok "單一樣本"; else warn "VCF 裡有 $ns 個樣本：後面的步驟假設只有孩子一個人"; fi
if bcftools view -h "$vcf" | grep -E '^##contig=<ID=chr' >/dev/null; then
  warn "染色體名稱有 chr 前綴：和 spec 的預期不同，參考序列也要用有 chr 的版本（先停下來確認）"
else
  ok "染色體名稱沒有 chr 前綴（bcftools 指令要寫 15:…，提交 CSV 才加 chr）"
fi
for t in PGT PID; do
  if bcftools view -h "$vcf" | grep -E "^##FORMAT=<ID=$t," >/dev/null; then ok "FORMAT 有 ${t}（讀段層級 phasing）"; else warn "FORMAT 沒有 ${t}：無法判斷 cis / trans，配對一律當未定相"; fi
done
info "提醒：INFO 的 AF 是這個人自己的等位基因比例（AN=2），不是族群頻率；族群頻率一律用 VEP 加上的 gnomAD"

say "HPO terms → $RESULTS/hpo_terms.tsv（孩子的表現型：放 results，不進 repo）"
set +e
python "$SCRIPTS/hpo_terms.py" extract --data "$DATA" --out "$RESULTS/hpo_terms.tsv" --obo "$HPO_OBO"
rc=$?
set -e
if [ "$rc" = 2 ]; then
  warn "資料裡沒有現成的 HPO 編號：手動把 HPO ID 和名稱填進 $RESULTS/hpo_terms.tsv（Tab 分隔），然後執行
    python track1/scripts/hpo_terms.py check --tsv $RESULTS/hpo_terms.tsv --obo $HPO_OBO"
elif [ "$rc" != 0 ]; then
  die "HPO 整理失敗（代碼 ${rc}）"
fi
say "完成。下一步：bash track1/scripts/04_qc.sh"
