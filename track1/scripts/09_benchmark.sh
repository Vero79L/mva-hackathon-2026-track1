#!/usr/bin/env bash
# 09 推廣性基準測試（Phase 2；事先登錄於 config/phase2_protocol.yml 的 benchmark 段）
# 只用公開資料：GIAB HG001–HG007 的基因體、PanelApp 基因組、ClinVar。完全不讀孩子的資料，
# 所以 Builder（Claude Code）也可以跑；中間檔放在 ~/mva-hackathon/public-benchmark（不是 work/）。
# 步驟：事先登錄檢查 → 下載 GIAB（約 0.9 GB）→ 每個疾病的基因範圍 → 背景基因體 → ClinVar 致病變異
#       → 共用目錄的 VEP 註解（和 05 相同參數）→ 模擬（凍結的 rank_panel.py）
# 輸出：track1/benchmark/output/benchmark.json、simulations.tsv（公開資料，跑完 commit）
# 用法：bash track1/scripts/09_benchmark.sh     約 1–2 小時（VEP 佔大半，視電腦而定）；完成的步驟會跳過，FORCE=1 全部重做
#       BENCH_WORKERS=6 bash track1/scripts/09_benchmark.sh   同時跑幾個 rank_panel.py（預設 4）
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
BENCH="${BENCH:-$MVA_ROOT/public-benchmark}"
mkdir -p "$BENCH/logs" "$BENCH/giab" "$BENCH/bg" "$BENCH/cat" "$BENCH/tmp"
LOGFILE="$BENCH/logs/09_benchmark_$(date +%Y%m%d-%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1
say "09_benchmark 開始於 $(date '+%Y-%m-%d %H:%M:%S')（紀錄：${LOGFILE}）"
conda_on mva
need bcftools "（先跑 01_setup_mac.sh）"
[ -s "$REF.fai" ] || die "缺參考序列：先跑 02_fetch_resources.sh fasta"
[ -s "$CLINVAR" ] || die "缺 ClinVar：先跑 02_fetch_resources.sh clinvar"
[ -s "$GTF" ] || die "缺 Ensembl GTF：先跑 02_fetch_resources.sh panel"
[ -f "$VEP_CACHE/homo_sapiens/${VEP_VERSION}_GRCh38/.mva_complete" ] || die "缺 VEP cache：先跑 02_fetch_resources.sh vep"
proto="$CONF/phase2_protocol.yml"
SAMPLES="HG001 HG002 HG003 HG004 HG005 HG006 HG007"

say "1/7 事先登錄檢查（規則、分析計畫、基準測試的程式和基因清單都要已 commit）"
check_preregistered
for f in config/phase2_protocol.yml benchmark/genelists_panelapp.tsv scripts/benchmark.py scripts/make_rules.py; do
  git -C "$T1" ls-files --error-unmatch "$f" >/dev/null 2>&1 || die "${f} 還沒 commit。先在終端機執行：
    cd \"$REPO\" && git add track1 && git commit -m \"Track 1: pre-register phase 2 analyses\""
  git -C "$T1" diff --quiet HEAD -- "$f" || die "${f} 有還沒 commit 的修改：先 commit（改計畫請把 version +1）"
done
pc="$(git -C "$CONF" log -1 --format='%h %ci' -- phase2_protocol.yml)"
ok "基準測試已事先登錄：commit ${pc}"

say "2/7 下載 GIAB 背景基因體（NISTv4.2.1 GRCh38，7 個，約 0.9 GB；已下載的會核對大小後跳過）"
GIAB="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release"
giab_dir() {
  case "$1" in
    HG001) echo "NA12878_HG001/NISTv4.2.1/GRCh38" ;;
    HG002) echo "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38" ;;
    HG003) echo "AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38" ;;
    HG004) echo "AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38" ;;
    HG005) echo "ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38" ;;
    HG006) echo "ChineseTrio/HG006_NA24694_father/NISTv4.2.1/GRCh38" ;;
    HG007) echo "ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/GRCh38" ;;
    *) die "不認識的 GIAB 樣本：$1" ;;
  esac
}
for s in $SAMPLES; do
  f="${s}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  fetch "$GIAB/$(giab_dir "$s")/$f" "$BENCH/giab/$f"
  fetch "$GIAB/$(giab_dir "$s")/$f.tbi" "$BENCH/giab/$f.tbi"
done

say "3/7 每個疾病的基因範圍（MVA 規則 + 7 個 PanelApp 基因組；基因 ±${PANEL_FLANK} bp）"
python "$SCRIPTS/benchmark.py" regions --gtf "$GTF" --genelists "$T1/benchmark/genelists_panelapp.tsv" \
  --meta "$T1/benchmark/genelists_panelapp.meta.json" --rules "$CONF/ranking_rules.yml" --protocol "$proto" \
  --flank "$PANEL_FLANK" --out "$BENCH/regions"

say "4/7 背景基因體：只取基因範圍、去掉 chr、和 05 相同的標準化"
for s in $SAMPLES; do
  in="$BENCH/giab/${s}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  out="$BENCH/bg/${s}.vcf.gz"
  if fresh "$BENCH/bg/${s}.gt.tsv" "$in" "$BENCH/regions/union.chr.bed" "$REF"; then ok "${s} 已完成"; continue; fi
  bcftools view -R "$BENCH/regions/union.chr.bed" "$in" -Ou \
    | bcftools annotate --rename-chrs "$BENCH/regions/chrmap.txt" -Ou - \
    | bcftools annotate -x INFO --set-id '%CHROM:%POS:%REF:%ALT' -Ou - \
    | bcftools norm -m -any -f "$REF" -c w -Ou - 2>> "$BENCH/logs/norm.log" \
    | bcftools view -e 'ALT="*"' -Oz -o "$out.part" -
  mv "$out.part" "$out"
  bcftools index -f -t "$out"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%FILTER\t%QUAL\t[%GT]\t[%AD]\t[%DP]\t[%GQ]\n' "$out" > "$BENCH/bg/${s}.gt.tsv.part"
  mv "$BENCH/bg/${s}.gt.tsv.part" "$BENCH/bg/${s}.gt.tsv"
  ok "${s}：$(wc -l < "$BENCH/bg/${s}.gt.tsv" | tr -d ' ') 個變異"
done

say "5/7 ClinVar 致病變異（P/LP、至少 1 星、在各疾病 tier A 隱性基因上、≤ 50 bp）"
spikes="$BENCH/cat/spikes.tsv"
if fresh "$spikes" "$CLINVAR" "$BENCH/regions/panels.json"; then ok "已完成"; else
  bcftools view -R "$BENCH/regions/union.bed" "$CLINVAR" -Ou \
    | bcftools view -i 'INFO/CLNSIG~"athogenic"' -Ou - \
    | bcftools norm -m -any -f "$REF" -c w -Ou - 2>> "$BENCH/logs/norm.log" \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/ORIGIN\n' \
    > "$BENCH/cat/spikes.raw.tsv"
  python "$SCRIPTS/benchmark.py" spikes --raw "$BENCH/cat/spikes.raw.tsv" --panels "$BENCH/regions/panels.json" \
    --out "$spikes.part" --vcf "$BENCH/cat/spikes.vcf" --fai "$REF.fai"
  bcftools sort -Oz -o "$BENCH/cat/spikes.vcf.gz" "$BENCH/cat/spikes.vcf" 2>/dev/null
  bcftools index -f -t "$BENCH/cat/spikes.vcf.gz"
  mv "$spikes.part" "$spikes"
fi

say "6/7 共用目錄：7 個背景與 ClinVar 變異的位點聯集 → VEP ${VEP_VERSION}（和 05 相同參數，--offline）"
sites="$BENCH/cat/catalog.sites.vcf.gz"
vepout="$BENCH/cat/catalog.vep.vcf.gz"
if fresh "$sites.tbi" "$BENCH/cat/spikes.vcf.gz" "$BENCH/bg/HG007.vcf.gz"; then ok "位點聯集已完成"; else
  lst="$BENCH/cat/sites.list"
  : > "$lst"
  for s in $SAMPLES; do
    bcftools view -G "$BENCH/bg/${s}.vcf.gz" -Oz -o "$BENCH/cat/sites_${s}.vcf.gz"
    bcftools index -f -t "$BENCH/cat/sites_${s}.vcf.gz"
    echo "$BENCH/cat/sites_${s}.vcf.gz" >> "$lst"
  done
  echo "$BENCH/cat/spikes.vcf.gz" >> "$lst"
  bcftools concat -a -D -f "$lst" -Ou | bcftools annotate -x INFO -Oz -o "$sites.part" -
  mv "$sites.part" "$sites"
  bcftools index -f -t "$sites"
fi
ok "$(bcftools index -n "$sites") 個不重複的變異"
if fresh "$vepout" "$sites" "$CLINVAR"; then ok "VEP 已完成"; else
  t0="$(date +%s)"
  conda_on vep
  vep --offline --cache --dir_cache "$VEP_CACHE" --cache_version "$VEP_VERSION" --assembly GRCh38 --fasta "$REF" \
      -i "$sites" -o "$vepout.part.vcf.gz" --vcf --compress_output bgzip \
      --everything --hgvsg --check_existing \
      --custom "file=$CLINVAR,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN" \
      --fork "$VEP_FORK" --force_overwrite --stats_text --no_progress
  conda_on mva
  mv "$vepout.part.vcf.gz" "$vepout"
  echo $(( ($(date +%s) - t0) / 60 )) > "$BENCH/cat/vep_minutes.txt"
fi
catalog="$BENCH/cat/catalog.tsv.gz"
if fresh "$catalog" "$vepout" "$BENCH/regions/panels.json"; then ok "可能合格的列已整理"; else
  fmt="$(python "$SCRIPTS/benchmark.py" format --splitvep "$CONF/splitvep_format.txt")"
  bcftools +split-vep "$vepout" -d -f "$fmt" \
    | python "$SCRIPTS/benchmark.py" catalog --panels "$BENCH/regions/panels.json" --rules "$CONF/ranking_rules.yml" \
        --out "$catalog.part.gz"
  mv "$catalog.part.gz" "$catalog"
fi

say "7/7 模擬：7 個背景 × 8 個疾病 × 15 個基因型 × 2 種模式（凍結的 rank_panel.py，同時 ${BENCH_WORKERS:-4} 個）"
python "$SCRIPTS/benchmark.py" simulate --regions "$BENCH/regions" --catalog "$catalog" --genotypes "$BENCH/bg" \
  --spikes "$spikes" --rules "$CONF/ranking_rules.yml" --protocol "$proto" --out "$T1/benchmark/output" \
  --workers "${BENCH_WORKERS:-4}" --tmp "$BENCH/tmp" --protocol-commit "${pc%% *}" --vep-version "$VEP_VERSION" \
  --clinvar-version "$(cat "$T1/env/clinvar_version.txt" 2>/dev/null || echo unknown)" \
  --annotation-minutes "$(cat "$BENCH/cat/vep_minutes.txt" 2>/dev/null || echo 0)"
( manifest 09_benchmark "protocol ${pc%% *}" ) 2>/dev/null || warn "無法寫入 run manifest（${RESULTS}）；結果不受影響"
say "完成。結果在 track1/benchmark/output/（公開資料）：
    cd \"$REPO\" && git add track1/benchmark && git commit -m \"Track 1: generalisation benchmark results\""
