#!/usr/bin/env bash
# 04 VCF 品質檢查與性別判斷（約 10 分鐘）
#   $WORK/qc/      完整統計（含逐染色體數字，比賽後刪）
#   $WORK/sex.txt  MALE / FEMALE / UNKNOWN（06 Exomiser 用；也可以在 config/local.env 設 SEX= 手動指定）
#   $T1/qc/qc_summary.txt  只有全基因體彙總數字，可以進 repo、寫進方法報告
# 用法：bash track1/scripts/04_qc.sh
. "$(dirname "$0")/lib.sh"
start_log 04_qc
conda_on mva
vcf="$(find_vcf)"
ensure_index "$vcf"
mkdir -p "$WORK/qc" "$T1/qc"

say "1/3 整體統計"
if fresh "$WORK/qc/stats.txt" "$vcf"; then ok "已完成（$WORK/qc/stats.txt）"; else
  bcftools stats -s - --threads "$THREADS" "$vcf" > "$WORK/qc/stats.txt.part"
  mv "$WORK/qc/stats.txt.part" "$WORK/qc/stats.txt"
fi
sn() { awk -F'\t' -v k="$1" '$1=="SN" && $3==k {print $4}' "$WORK/qc/stats.txt"; }   # SN [2]id [3]key [4]value
records="$(sn 'number of records:')"
snps="$(sn 'number of SNPs:')"
indels="$(sn 'number of indels:')"
multi="$(sn 'number of multiallelic sites:')"
tstv="$(awk -F'\t' '$1=="TSTV" {print $5; exit}' "$WORK/qc/stats.txt")"
# PSC：[4] nRefHom [5] nNonRefHom [6] nHets [10] average depth
read -r homalt hets depth <<EOF
$(awk -F'\t' '$1=="PSC" {print $5, $6, $10; exit}' "$WORK/qc/stats.txt")
EOF
hethom="$(awk -v h="$hets" -v o="$homalt" 'BEGIN { if (o > 0) printf "%.2f", h / o; else print "NA" }')"
bcftools query -f '%FILTER\n' "$vcf" | sort | uniq -c | sort -nr > "$WORK/qc/filter_counts.txt"
bcftools query -f '[%GT]\n' "$vcf" | sort | uniq -c | sort -nr > "$WORK/qc/gt_counts.txt"
bcftools query -f '%CHROM\n' "$vcf" | uniq -c > "$WORK/qc/per_chrom.txt"
pass="$(awk '$2=="PASS" {print $1}' "$WORK/qc/filter_counts.txt")"
passpct="$(awk -v p="${pass:-0}" -v t="$records" 'BEGIN { if (t > 0) printf "%.1f", 100 * p / t; else print "NA" }')"

say "2/3 性別判斷（X 非 PAR 區的 het/hom-alt 比、Y 的 call 數）"
xr="X:2781480-155701382"
contigs="$(bcftools view -h "$vcf" | sed -nE 's/^##contig=<ID=([^,>]+).*/\1/p')"
if printf '%s\n' "$contigs" | grep -x 'chrX' >/dev/null; then xr="chr$xr"; yc=chrY; else yc=Y; fi
xhet="$(bcftools view -H -r "$xr" -i 'GT="het"' "$vcf" | wc -l | tr -d ' ')"
xhom="$(bcftools view -H -r "$xr" -i 'GT="AA"' "$vcf" | wc -l | tr -d ' ')"
ycalls="$( (bcftools view -H -r "$yc" -i 'GT="alt"' "$vcf" 2>/dev/null || true) | wc -l | tr -d ' ')"
xratio="$(awk -v h="$xhet" -v o="$xhom" 'BEGIN { if (o > 0) printf "%.2f", h / o; else print "NA" }')"
sex=UNKNOWN
if [ "$xratio" != NA ]; then
  sex="$(awk -v r="$xratio" -v y="$ycalls" 'BEGIN { if (r < 0.5 && y >= 1000) print "MALE"; else if (r >= 0.8 && y < 500) print "FEMALE"; else print "UNKNOWN" }')"
fi
if [ -n "${SEX:-}" ]; then
  info "config/local.env 指定 SEX=${SEX}（自動判斷結果：${sex}）"
  sex="$SEX"
fi
echo "$sex" > "$WORK/sex.txt"
echo "chrX nonPAR het=$xhet homalt=$xhom ratio=$xratio; chrY alt calls=$ycalls; inferred=$sex" > "$WORK/qc/sex_check.txt"
info "$(cat "$WORK/qc/sex_check.txt")"
if [ "$sex" = UNKNOWN ]; then
  warn "判斷不出來：對照資料集說明裡的性別，在 config/local.env 加一行 SEX=MALE 或 SEX=FEMALE"
else
  ok "性別：${sex}（男性 X 的比值通常 < 0.3 且 Y 有上千個 call；請和資料集說明對照，不一致就停下來查）"
fi

say "3/3 摘要 → $T1/qc/qc_summary.txt（只有全基因體彙總數字）"
cat > "$T1/qc/qc_summary.txt" <<EOF
# Track 1 VCF QC summary (genome-wide aggregates only; generated $(date '+%Y-%m-%d'))
records	$records
snps	$snps
indels	$indels
multiallelic_sites	$multi
ts_tv	$tstv
het_hom_ratio	$hethom
mean_depth	$depth
pass_percent	$passpct
EOF
cat "$T1/qc/qc_summary.txt"
echo
awk -v r="$records" -v t="$tstv" -v h="$hethom" -v p="$passpct" 'BEGIN {
  printf "  判讀（30x WGS 參考範圍）\n"
  printf "    總變異數 %s（400–550 萬）%s\n", r, (r >= 4000000 && r <= 5500000) ? "" : "  ← 偏離，確認是否為完整 WGS"
  printf "    Ti/Tv %s（2.0–2.1）%s\n", t, (t >= 1.95 && t <= 2.15) ? "" : "  ← 偏離：假陽性或過濾問題"
  printf "    het/hom %s（1.4–2.0）%s\n", h, (h == "NA" || (h >= 1.4 && h <= 2.0)) ? "" : "  ← 偏離：偏低可能有大片段同型合子（近親背景）"
  printf "    PASS %s%%\n", p
}'
say "完成。下一步：bash track1/scripts/05_annotate.sh"
