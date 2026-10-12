#!/usr/bin/env bash
# 02 下載公開參考資源到 ${RES}（比賽結束後可以保留）。斷線重跑會續傳；已完成的步驟自動跳過。
# 用法：bash track1/scripts/02_fetch_resources.sh [步驟...]
#   步驟：scorer hpo fasta clinvar panel exomiser vep（不給 = 全部，依這個順序；大檔放最後）
# 空間：全部約 150 GB（VEP cache 解壓尖峰約 60 GB、Exomiser 資料解壓尖峰約 2 倍）。建議晚上接電源掛著跑。
. "$(dirname "$0")/lib.sh"
stay_awake "$@"
start_log 02_fetch_resources
conda_on mva
need bgzip
need samtools
mkdir -p "$RES"

URL_FASTA="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta.gz"
URL_VEP="https://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/indexed_vep_cache/homo_sapiens_vep_${VEP_VERSION}_GRCh38.tar.gz"
URL_GTF="https://ftp.ensembl.org/pub/release-${VEP_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${VEP_VERSION}.gtf.gz"
URL_CLINVAR="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
URL_EXOMISER="https://github.com/exomiser/Exomiser/releases/download/$EXOMISER_VERSION/exomiser-cli-$EXOMISER_VERSION-distribution.zip"
URL_EXOMISER_DATA="https://g-879a9f.f5dc97.75bc.dn.glob.us/data"
URL_HPO="https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo"
URL_SCORER="https://huggingface.co/spaces/$HF_SPACE/raw/main/evaluation.py"

check_contigs() {   # 參考序列的主要染色體名稱與長度要和 VCF 檔頭完全一致
  local vcf
  vcf="$(find_vcf 2>/dev/null)" || { info "（還沒有 VCF，等 03_fetch_data.sh 之後再比對染色體）"; return 0; }
  local tab
  tab="$(printf '\t')"
  bcftools view -h "$vcf" | sed -nE "s/^##contig=<ID=([^,>]+),length=([0-9]+).*/\\1${tab}\\2/p" \
    | grep -E "^(chr)?([0-9]+|X|Y|M|MT)${tab}" | sort > "$RES/.contigs_vcf.tsv" || true
  cut -f1,2 "$REF.fai" | grep -E "^(chr)?([0-9]+|X|Y|M|MT)${tab}" | sort > "$RES/.contigs_ref.tsv" || true
  if [ ! -s "$RES/.contigs_vcf.tsv" ]; then warn "VCF 檔頭沒有 contig 長度，略過比對"; return 0; fi
  if diff "$RES/.contigs_vcf.tsv" "$RES/.contigs_ref.tsv" >/dev/null; then
    ok "參考序列與 VCF 的 $(wc -l < "$RES/.contigs_vcf.tsv" | tr -d ' ') 條主要染色體名稱、長度完全一致"
  else
    diff "$RES/.contigs_vcf.tsv" "$RES/.contigs_ref.tsv" | head -10
    die "參考序列和 VCF 的染色體名稱或長度不一致：不要往下做，先確認 VCF 檔頭的 ##reference"
  fi
}

get_fasta() {
  say "GRCh38 參考序列（no-alt analysis set，去掉 chr 前綴）"
  if [ -s "$REF.fai" ] && [ -s "$REF.gzi" ] && [ -s "$REF_PLAIN.fai" ]; then ok "已完成"; check_contigs; return 0; fi
  need_space "$RES" 12
  local gz
  gz="$RES/$(basename "$URL_FASTA")"
  fetch "$URL_FASTA" "$gz"
  info "轉成無 chr 前綴的 bgzip 檔並建索引（約 5–10 分鐘）"
  gzip -dc "$gz" | sed -E 's/^>chr/>/' | bgzip -@ "$THREADS" > "$REF.part"
  mv "$REF.part" "$REF"
  samtools faidx "$REF"
  bgzip -dc "$REF" > "$REF_PLAIN.part"
  mv "$REF_PLAIN.part" "$REF_PLAIN"
  samtools faidx "$REF_PLAIN"
  rm -f "$gz"
  ok "$(cut -f1 "$REF.fai" | head -25 | tr '\n' ' ')"
  check_contigs
}

get_clinvar() {
  say "ClinVar（最新版，比 VEP cache 內建的新）"
  if [ -s "$CLINVAR" ] && [ "${FORCE:-0}" != 1 ]; then
    ok "已存在：$(bcftools view -h "$CLINVAR" | grep -m1 '^##fileDate' || true)（要更新就 FORCE=1 重跑）"
  else
    rm -f "$CLINVAR" "$CLINVAR.tbi"
    fetch "$URL_CLINVAR" "$CLINVAR"
  fi
  [ -s "$CLINVAR.tbi" ] || bcftools index -t "$CLINVAR"
  mkdir -p "$T1/env"
  bcftools view -h "$CLINVAR" | grep -m1 '^##fileDate' > "$T1/env/clinvar_version.txt" || true
  ok "ClinVar $(cat "$T1/env/clinvar_version.txt")"
}

get_panel() {
  say "MVA 基因組 BED（13 個基因 ±${PANEL_FLANK} bp；基因名單來自 ranking_rules.yml）"
  fetch "$URL_GTF" "$GTF"
  local re aliases n got
  re="$(python "$SCRIPTS/panel_regions.py" genes --regex)"
  aliases="$(python "$SCRIPTS/panel_regions.py" genes --aliases)"
  n="$(python "$SCRIPTS/panel_regions.py" genes --count)"
  gzip -dc "$GTF" | awk -F'\t' -v OFS='\t' -v flank="$PANEL_FLANK" -v re="gene_name \"($re)\";" -v al="$aliases" '
    BEGIN { k = split(al, a, ","); for (i = 1; i <= k; i++) { split(a[i], p, "="); map[p[1]] = p[2] } }
    $3 == "gene" && $9 ~ re {
      match($9, /gene_name "[^"]+"/); g = substr($9, RSTART + 11, RLENGTH - 12)
      if (g in map) g = map[g]
      s = $4 - 1 - flank; if (s < 0) s = 0
      print $1, s, $5 + flank, g
    }' | sort -k1,1 -k2,2n > "$PANEL_BED.part"
  mv "$PANEL_BED.part" "$PANEL_BED"
  got="$(cut -f4 "$PANEL_BED" | sort -u | wc -l | tr -d ' ')"
  [ "$got" = "$n" ] || die "BED 裡只有 $got 個基因，應該是 $n 個：$(cut -f4 "$PANEL_BED" | sort -u | tr '\n' ' ')。
  少掉的基因多半是 Ensembl 用了別的名稱，把舊名加進 ranking_rules.yml 的 genes.aliases 再跑一次。"
  ok "$got 個基因：$(cut -f4 "$PANEL_BED" | sort -u | tr '\n' ' ')"
}

get_hpo() {
  say "HPO 本體 hp.obo（檢查 HPO 編號是否過時、補上名稱）"
  if [ "${FORCE:-0}" = 1 ]; then rm -f "$HPO_OBO"; fi
  fetch "$URL_HPO" "$HPO_OBO"
  ok "$(grep -m1 '^data-version' "$HPO_OBO" || true)"
}

get_scorer() {
  say "官方評分程式 evaluation.py（08_submission.sh 自我評分用）"
  rm -f "$SCORER"          # 每次都抓最新版，主辦方可能更新
  fetch "$URL_SCORER" "$SCORER"
  grep -E '^def (load_submission|score_proband)' "$SCORER" >/dev/null \
    || warn "evaluation.py 裡沒看到 load_submission / score_proband：官方可能改了介面，08 會改用內建檢查並提示"
  ok "sha256 $(sha256 "$SCORER" | cut -c1-16)…"
}

get_vep() {
  say "VEP $VEP_VERSION cache（GRCh38；下載 25–30 GB，解壓尖峰約 60 GB）"
  local d="$VEP_CACHE/homo_sapiens/${VEP_VERSION}_GRCh38" tgz
  if [ -f "$d/.mva_complete" ]; then ok "已完成"; return 0; fi
  tgz="$VEP_CACHE/$(basename "$URL_VEP")"
  [ -s "$tgz" ] || need_space "$RES" 65
  fetch "$URL_VEP" "$tgz"
  info "解壓（約 30–60 分鐘）"
  tar -xzf "$tgz" -C "$VEP_CACHE"
  [ -d "$d" ] || die "解壓後找不到 $d"
  touch "$d/.mva_complete"
  rm -f "$tgz"
  ok "VEP cache：$d"
}

get_exomiser() {
  say "Exomiser $EXOMISER_VERSION + 資料 ${EXOMISER_DATA_VERSION}（hg38 + phenotype）"
  local base="$RES/exomiser" jar="$EXOMISER_DIR/exomiser-cli-$EXOMISER_VERSION.jar" zip part name
  mkdir -p "$base" "$EXOMISER_DATA"
  if [ ! -s "$jar" ]; then
    zip="$base/exomiser-cli-$EXOMISER_VERSION-distribution.zip"
    fetch "$URL_EXOMISER" "$zip"
    unzip_big "$zip" "$base"
    [ -s "$jar" ] || die "解壓後找不到 ${jar}（GitHub 下載失敗時，到 Exomiser 的 releases 頁面手動下載 distribution zip 放到 ${base}）"
    rm -f "$zip"
  fi
  ok "$(basename "$jar")"
  for part in hg38 phenotype; do
    name="${EXOMISER_DATA_VERSION}_$part"
    if [ -f "$EXOMISER_DATA/$name/.mva_complete" ]; then ok "$name 已完成"; continue; fi
    zip="$base/$name.zip"
    if [ ! -s "$zip" ]; then
      if [ "$part" = hg38 ]; then need_space "$EXOMISER_DATA" 80; else need_space "$EXOMISER_DATA" 20; fi
    fi
    fetch "$URL_EXOMISER_DATA/$name.zip" "$zip"
    info "解壓 ${name}（可能要 20–60 分鐘）"
    unzip_big "$zip" "$EXOMISER_DATA"
    [ -d "$EXOMISER_DATA/$name" ] || die "解壓後找不到 $EXOMISER_DATA/$name"
    touch "$EXOMISER_DATA/$name/.mva_complete"
    rm -f "$zip"
  done
  [ -f "$EXOMISER_DIR/application.properties.bak" ] || cp "$EXOMISER_DIR/application.properties" "$EXOMISER_DIR/application.properties.bak" 2>/dev/null || true
  cat > "$EXOMISER_DIR/application.properties" <<EOF
exomiser.data-directory=$EXOMISER_DATA
exomiser.hg38.data-version=$EXOMISER_DATA_VERSION
exomiser.phenotype.data-version=$EXOMISER_DATA_VERSION
EOF
  ok "application.properties 已指向 ${EXOMISER_DATA}（原檔備份為 application.properties.bak）"
  need java "（Java 21：先跑 01_setup_mac.sh）"
  (cd "$EXOMISER_DIR" && java -jar "$jar" --help 2>&1 | head -5) || true
  printf 'exomiser-cli %s, data %s\n' "$EXOMISER_VERSION" "$EXOMISER_DATA_VERSION" > "$T1/env/exomiser.txt"
}

steps="${*:-scorer hpo fasta clinvar panel exomiser vep}"
for s in $steps; do
  case "$s" in
    fasta) get_fasta ;;
    clinvar) get_clinvar ;;
    panel) get_panel ;;
    hpo) get_hpo ;;
    scorer) get_scorer ;;
    vep) get_vep ;;
    exomiser) get_exomiser ;;
    *) die "不認識的步驟：${s}（可用：scorer hpo fasta clinvar panel exomiser vep）" ;;
  esac
done
say "完成：${steps}。剩餘空間 $(free_gb "$RES") GB。下一步：bash track1/scripts/03_fetch_data.sh"
