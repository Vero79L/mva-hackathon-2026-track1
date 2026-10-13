#!/usr/bin/env bash
# 06 Exomiser（路線 A）：不預設任何基因，用孩子的 HPO 表現型對全基因體變異打分排序。跑兩次：
#   exome    官方預設 preset（評估編碼區與剪接區）
#   introns  15.x 新增的含內含子 preset（examples/preset-exome-analysis-with-introns.yml，有這個檔才跑）
# 輸出含逐位點資料，放 $WORK/exomiser（比賽後刪）；基因層級摘要放 $RESULTS/exomiser_summary.md。
# 每次約 20–60 分鐘；記憶體不足（OutOfMemoryError）就在 config/local.env 把 JAVA_XMX 調低。
# 用法：bash track1/scripts/06_exomiser.sh
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
start_log 06_exomiser
conda_on mva
need java "（Java 21：先跑 01_setup_mac.sh）"
jar="$EXOMISER_DIR/exomiser-cli-$EXOMISER_VERSION.jar"
[ -s "$jar" ] || die "缺 Exomiser：先跑 02_fetch_resources.sh exomiser"
for p in hg38 phenotype; do
  [ -f "$EXOMISER_DATA/${EXOMISER_DATA_VERSION}_$p/.mva_complete" ] || die "缺 Exomiser 資料 ${EXOMISER_DATA_VERSION}_${p}：先跑 02_fetch_resources.sh exomiser"
done
norm="$WORK/norm.vcf.gz"
[ -s "$norm" ] || die "缺標準化 VCF：先跑 05_annotate.sh"
[ -s "$RESULTS/hpo_terms.tsv" ] || die "缺 HPO 清單：先跑 03_fetch_data.sh（或手動建立 $RESULTS/hpo_terms.tsv）"
if ! (check_preregistered) >/dev/null 2>&1; then
  warn "排序規則（ranking_rules.yml）還沒 commit。Exomiser 的結果也算「結果」：建議先 commit 規則再往下看（見 README）"
fi
sex="$(cat "$WORK/sex.txt" 2>/dev/null || echo UNKNOWN)"     # 04_qc.sh 已套用 sex_override.txt 或 SEX
ex="$WORK/exomiser"
mkdir -p "$ex"

say "1/3 phenopacket（從 HPO 清單產生；subject.id 自動取 VCF 樣本名稱）"
pp="$ex/proband.phenopacket.yml"
python "$SCRIPTS/make_phenopacket.py" --hpo "$RESULTS/hpo_terms.tsv" --vcf "$norm" --sex "$sex" \
  --obo "$HPO_OBO" --id "$PROBAND" --out "$pp.new"
if [ -s "$pp" ] && cmp -s <(grep -v '^  created:' "$pp") <(grep -v '^  created:' "$pp.new"); then
  rm -f "$pp.new"; ok "phenopacket 沒有變（保留原檔，避免不必要的重跑）"
else
  mv "$pp.new" "$pp"
fi

say "2/3 Exomiser ${EXOMISER_VERSION}（資料 ${EXOMISER_DATA_VERSION}，-Xmx${JAVA_XMX}）"
cd "$EXOMISER_DIR"
sub=""
if java -jar "$jar" --help 2>&1 | grep -E '(^|[[:space:]])analyse([[:space:]]|$)' >/dev/null; then sub="analyse"; fi

run_exomiser() {   # run_exomiser <名稱> <分析參數...>
  local name="$1" out
  shift
  out="$ex/${PROBAND}_$name.genes.tsv"
  if fresh "$out" "$norm" "$pp" "$EXOMISER_DIR/application.properties"; then ok "$name 已完成（${out}）"; return 0; fi
  cat > "$ex/output-$name.yml" <<EOF
outputContributingVariantsOnly: false
numGenes: 0
outputDirectory: "$ex"
outputFileName: "${PROBAND}_$name"
outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT]
EOF
  info "${name}：開始於 $(date '+%H:%M')"
  # shellcheck disable=SC2086
  java "-Xmx$JAVA_XMX" -jar "$jar" $sub "$@" --sample "$pp" --vcf "$norm" --assembly hg38 --output "$ex/output-$name.yml"
  [ -s "$out" ] || die "Exomiser $name 沒有產生 ${out}：看上面的錯誤訊息（記憶體不足就調低 JAVA_XMX）"
  ok "$name 完成於 $(date '+%H:%M')"
}

run_exomiser exome --preset exome
if [ -f examples/preset-exome-analysis-with-introns.yml ]; then
  run_exomiser introns --analysis examples/preset-exome-analysis-with-introns.yml
else
  warn "這個版本沒有 examples/preset-exome-analysis-with-introns.yml：只跑 exome preset"
fi
cd - >/dev/null

say "3/3 基因層級摘要 → $RESULTS/exomiser_summary.md"
python "$SCRIPTS/exomiser_summary.py" --dir "$ex" --out "$RESULTS/exomiser_summary.md" --rules "$CONF/ranking_rules.yml"
info "完整報告（含變異）：open $ex/${PROBAND}_exome.html"
say "完成。下一步：bash track1/scripts/07_rank.sh"
