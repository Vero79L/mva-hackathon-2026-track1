#!/usr/bin/env bash
# 05 標準化 → 抽出 MVA 基因組 → VEP 完整註解（完全離線）→ 表格 → SpliceAI（選配）→ SV / CNV 比對
# 輸出都在 ${WORK}（含基因型，比賽後刪）：norm.vcf.gz、panel.vcf.gz、panel.vep.vcf.gz、panel.vep.tsv、
# panel.spliceai.tsv、panel_sv.tsv。已完成而且輸入沒變的步驟會跳過；FORCE=1 全部重做。
# 用法：bash track1/scripts/05_annotate.sh          約 20–40 分鐘（SpliceAI 佔大半）
#       SKIP_SPLICEAI=1 bash track1/scripts/05_annotate.sh
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
start_log 05_annotate
conda_on mva
vcf="$(find_vcf)"
ensure_index "$vcf"
[ -s "$REF.fai" ] || die "缺參考序列：先跑 02_fetch_resources.sh fasta"
[ -s "$PANEL_BED" ] || die "缺 MVA 基因組 BED：先跑 02_fetch_resources.sh panel"
[ -s "$CLINVAR" ] || die "缺 ClinVar：先跑 02_fetch_resources.sh clinvar"
[ -f "$VEP_CACHE/homo_sapiens/${VEP_VERSION}_GRCh38/.mva_complete" ] || die "缺 VEP cache：先跑 02_fetch_resources.sh vep"
mkdir -p "$WORK"
RULES="$CONF/ranking_rules.yml"

say "1/6 標準化（多等位基因拆開、indel 左對齊；評分是逐字比對 POS/REF/ALT，這一步直接影響分數）"
norm="$WORK/norm.vcf.gz"
if fresh "$norm.tbi" "$vcf" "$REF"; then ok "已完成（${norm}）"; else
  regions="$(main_contigs "$vcf")"
  [ -n "$regions" ] || die "VCF 和參考序列沒有共同的主要染色體名稱：檢查 02 的參考序列"
  if ! bcftools view -r "$regions" "$vcf" -Ou \
      | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou - \
      | bcftools norm -m -any -f "$REF" -c w -Ou - 2> "$WORK/norm.log" \
      | bcftools view -e 'ALT="*"' -Oz -o "$norm.part" --threads "$THREADS" -; then
    tail -5 "$WORK/norm.log" >&2
    die "標準化失敗：訊息見上面（完整內容在 $WORK/norm.log）"
  fi
  mv "$norm.part" "$norm"
  bcftools index -f -t "$norm"
fi
info "$(grep -E 'total/split' "$WORK/norm.log" | tail -1 || true)"
mism="$(grep -ci 'mismatch' "$WORK/norm.log" || true)"
if [ "${mism:-0}" -gt 100 ]; then
  die "REF 和參考序列不符的有 $mism 筆：參考序列版本不對，回 02_fetch_resources.sh fasta 檢查（不要往下做）"
fi
ok "REF 不符 ${mism:-0} 筆（應接近 0）；原始寫法存在 ID 欄，例如 15:123:AT:A,ATT"

say "2/6 抽出 MVA 基因組（13 個基因 ±${PANEL_FLANK} bp + ranking_rules.yml 的已知等位基因；含非 PASS、含非編碼區）"
python "$SCRIPTS/panel_regions.py" bed --rules "$RULES" --bed "$PANEL_BED" --out "$WORK/panel_regions.bed.new"
if cmp -s "$WORK/panel_regions.bed.new" "$WORK/panel_regions.bed"; then
  rm -f "$WORK/panel_regions.bed.new"      # 範圍沒變：保留舊檔，後面的步驟才不會白白重跑
else
  mv "$WORK/panel_regions.bed.new" "$WORK/panel_regions.bed"
fi
panel="$WORK/panel.vcf.gz"
if fresh "$panel.tbi" "$norm" "$WORK/panel_regions.bed"; then ok "已完成"; else
  bcftools view -R "$WORK/panel_regions.bed" "$norm" -Oz -o "$panel.part"
  mv "$panel.part" "$panel"
  bcftools index -f -t "$panel"
fi
ok "$(bcftools view -H "$panel" | wc -l | tr -d ' ') 個變異"

say "3/6 VEP $VEP_VERSION 完整註解（--offline：座標不會離開這台電腦）"
vepout="$WORK/panel.vep.vcf.gz"
if fresh "$vepout" "$panel" "$CLINVAR"; then ok "已完成"; else
  conda_on vep
  vep --offline --cache --dir_cache "$VEP_CACHE" --cache_version "$VEP_VERSION" --assembly GRCh38 --fasta "$REF" \
      -i "$panel" -o "$vepout.part.vcf.gz" --vcf --compress_output bgzip \
      --everything --hgvsg --check_existing \
      --custom "file=$CLINVAR,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN" \
      --fork "$VEP_FORK" --force_overwrite --stats_text --no_progress
  mv "$vepout.part.vcf.gz" "$vepout"
  [ -f "$vepout.part.vcf.gz_summary.txt" ] && mv "$vepout.part.vcf.gz_summary.txt" "$WORK/panel.vep_summary.txt"
  conda_on mva
fi

say "4/6 轉成表格（每個轉錄本一行；排序時再挑 MANE / canonical）"
tsv="$WORK/panel.vep.tsv"
fmt="$(cat "$CONF/splitvep_format.txt")"
avail="$(bcftools +split-vep "$vepout" -l | awk '{print $2}')"
missing=""
for f in $(printf '%s' "$fmt" | grep -oE '%[A-Za-z_]+' | sed 's/^%//'); do
  case "$f" in CHROM|POS|ID|REF|ALT|FILTER|QUAL|GT|AD|DP|GQ|PGT|PID) continue ;; esac
  printf '%s\n' "$avail" | grep -x "$f" >/dev/null || missing="$missing $f"
done
[ -z "$missing" ] || die "VEP 輸出裡沒有這些欄位：$missing
  用  bcftools +split-vep $vepout -l  看實際名稱，改 config/splitvep_format.txt（rank_panel.py 的 COLS 也要對應）"
if fresh "$tsv" "$vepout" "$CONF/splitvep_format.txt"; then ok "已完成"; else
  bcftools +split-vep "$vepout" -d -f "$fmt" > "$tsv.part"
  mv "$tsv.part" "$tsv"
fi
ok "$(wc -l < "$tsv" | tr -d ' ') 行（變異 × 轉錄本）"

say "5/6 SpliceAI（評估內含子、UTR、同義變異會不會破壞剪接；-D ${SPLICEAI_DIST}）"
sai="$WORK/panel.spliceai.tsv"
if [ "${SKIP_SPLICEAI:-0}" = 1 ]; then
  info "SKIP_SPLICEAI=1：略過（排序時 SpliceAI 分數當 0）"
elif ! has_env spliceai; then
  warn "沒有 spliceai 環境：略過（排序時 SpliceAI 分數當 0；第二階段補）"
elif fresh "$sai" "$panel"; then
  ok "已完成"
else
  [ -s "$REF_PLAIN.fai" ] || die "缺未壓縮的參考序列：重跑 02_fetch_resources.sh fasta"
  conda_on spliceai
  spliceai -I "$panel" -O "$WORK/panel.spliceai.vcf" -R "$REF_PLAIN" -A grch38 -D "$SPLICEAI_DIST"
  conda_on mva
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SpliceAI\n' "$WORK/panel.spliceai.vcf" > "$sai.part"
  mv "$sai.part" "$sai"
  ok "$(awk -F'\t' '$5 != "."' "$sai" | wc -l | tr -d ' ') 個變異有 SpliceAI 分數"
fi

say "6/6 SV / CNV（在 MVA 基因上看過一遍；跨進基因的大片段也算）"
svout="$WORK/panel_sv.tsv"
: > "$svout"
sv="$(sv_vcfs)"
if [ -z "$sv" ]; then
  info "資料裡沒有 SV / CNV VCF"
else
  printf '%s\n' "$sv" | while IFS= read -r f; do
    ensure_index "$f" 2>/dev/null || { warn "無法建立索引：$f"; continue; }
    bed="$WORK/panel_regions.bed"
    if bcftools view -h "$f" | grep -E '^##contig=<ID=chr' >/dev/null; then
      sed 's/^/chr/' "$WORK/panel_regions.bed" > "$WORK/panel_regions.chr.bed"
      bed="$WORK/panel_regions.chr.bed"
    fi
    bcftools view --regions-overlap 1 -R "$bed" "$f" 2>/dev/null \
      | bcftools query -f "$(basename "$f")\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO\t[%GT]\n" >> "$svout" || true
  done
  n="$(wc -l < "$svout" | tr -d ' ')"
  if [ "$n" -gt 0 ]; then
    warn "$n 個 SV / CNV 和 MVA 基因組重疊：一定要人工看過（第二個等位基因可能就是它）→ $svout"
  else
    ok "沒有 SV / CNV 落在 MVA 基因組"
  fi
fi
say "完成。下一步：bash track1/scripts/06_exomiser.sh（排序規則要先 commit，見 README）"
