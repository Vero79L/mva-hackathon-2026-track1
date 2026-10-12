#!/usr/bin/env bash
# 00 檢查：Mac 規格、資料安全設定、工具、資源與進度。只讀不改，隨時可以跑；
# 唯一的例外是第一次執行時，依記憶體大小建立 config/local.env（個人設定，不進 git）。
# 用法：bash track1/scripts/00_check.sh
. "$(dirname "$0")/lib.sh"
set +e    # 檢查腳本：某一項沒過也要繼續往下看

chk() {   # chk <說明> <指令...>：指令成功打 ✓，失敗打 !
  local label="$1"
  shift
  if "$@"; then ok "$label"; else warn "$label"; fi
}
filevault_on() { fdesetup status 2>/dev/null | grep "On" >/dev/null; }
tm_excluded()  { tmutil isexcluded "$1" 2>/dev/null | grep Excluded >/dev/null; }
has_brew()     { command -v brew >/dev/null 2>&1 || [ -x /opt/homebrew/bin/brew ] || [ -x /usr/local/bin/brew ]; }
has_java21()   { if [ -x /usr/libexec/java_home ]; then /usr/libexec/java_home -v 21+ >/dev/null 2>&1; else command -v java >/dev/null 2>&1; fi; }
rules_done()   { [ -n "$(git -C "$CONF" log -1 --format=%h -- ranking_rules.yml 2>/dev/null)" ]; }

say "這台電腦"
arch="$(uname -m)"
if is_mac; then
  os="macOS $(sw_vers -productVersion 2>/dev/null)"
  mem_gb=$(( $(sysctl -n hw.memsize) / 1073741824 ))
  cores="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)"
else
  os="$(uname -s) $(uname -r)"
  mem_gb=$(( $(awk '/MemTotal/ {print $2}' /proc/meminfo) / 1048576 ))
  cores="$(getconf _NPROCESSORS_ONLN)"
fi
info "$os · $arch · ${mem_gb} GB 記憶體 · ${cores} 個效能核心"
if [ "$arch" = arm64 ]; then info "Apple Silicon：個別套件沒有原生版本時，01_setup_mac.sh 會自動改用 Rosetta（osx-64）"; fi

if [ "$mem_gb" -lt 12 ]; then xmx=6g; fork=2
elif [ "$mem_gb" -lt 24 ]; then xmx=10g; fork=4
else xmx=16g; fork=$(( cores > 8 ? 8 : cores )); fi
if [ ! -f "$CONF/local.env" ]; then
  cat > "$CONF/local.env" <<EOF
# shellcheck shell=bash
# 個人設定（不進 git）。依這台 Mac 的規格自動產生：${mem_gb} GB 記憶體、${cores} 個效能核心。
# 寫成 : "\${變數:=值}"：執行前 export 的同名變數仍然優先。
: "\${JAVA_XMX:=$xmx}"
: "\${VEP_FORK:=$fork}"
: "\${THREADS:=$fork}"
# 其他可覆蓋的設定（見 config/paths.env），需要時把開頭的 # 拿掉：
# : "\${RES:=/Volumes/MVA-SSD/resources}"      # 參考資源放外接 SSD
# : "\${VCF:=$DATA/<主 VCF 檔名>.vcf.gz}"       # 資料夾裡有不只一個 VCF 時指定
# : "\${SEX:=MALE}"                             # 04_qc.sh 判斷不出性別時手動指定
EOF
  ok "已建立 config/local.env：JAVA_XMX=${xmx}、VEP_FORK=$fork"
else
  info "config/local.env 已存在（JAVA_XMX=${JAVA_XMX}、VEP_FORK=${VEP_FORK}）"
fi

say "資料安全"
if is_mac; then chk "FileVault：$(fdesetup status 2>/dev/null | head -1)" filevault_on; fi
case "$MVA_ROOT" in
  "$HOME/Desktop"*|"$HOME/Documents"*|"$HOME/Library/Mobile Documents"*)
    warn "$MVA_ROOT 在桌面 / 文件 / iCloud 底下：可能被同步到 iCloud，請移到家目錄底下" ;;
  *) ok "$MVA_ROOT 不在 iCloud 同步範圍" ;;
esac
if command -v tmutil >/dev/null 2>&1; then
  for d in "$DATA" "$MVA_ROOT/work" "$WORK"; do
    if [ -d "$d" ]; then chk "Time Machine 不備份 $d" tm_excluded "$d"; fi
  done
fi
chk "git pre-commit 防護（01_setup_mac.sh 會裝）" test -x "$REPO/.git/hooks/pre-commit"
chk "track1/.gitignore" test -f "$T1/.gitignore"

say "工具"
if is_mac; then chk "Xcode 命令列工具" xcode-select -p; fi
chk "Homebrew" has_brew
chk "Java 21 以上（Exomiser 需要）" has_java21
if b="$(conda_base)"; then
  ok "conda：$b"
  for e in mva vep spliceai; do
    if has_env "$e"; then ok "conda 環境 $e"
    elif [ "$e" = spliceai ]; then info "conda 環境 spliceai 沒有（選配）"
    else warn "conda 環境 $e 沒有"; fi
  done
else
  warn "沒有 conda（01_setup_mac.sh 會裝 Miniforge）"
fi

say "比賽資料（${DATA}）"
if [ -d "$DATA" ]; then
  if v="$(find_vcf 2>/dev/null)"; then ok "主 VCF：${v#"$DATA"/}"; else warn "還沒有主 VCF（03_fetch_data.sh）"; fi
  chk "HPO 清單 results/track1/hpo_terms.tsv（03）" test -s "$RESULTS/hpo_terms.tsv"
else
  warn "$DATA 不存在"
fi

say "參考資源（${RES}，所在磁碟剩 $(free_gb "$RES") GB；第一階段建議至少 200 GB）"
chk "GRCh38 參考序列（02 fasta）" test -s "$REF.fai"
chk "VEP $VEP_VERSION cache（02 vep）" test -f "$VEP_CACHE/homo_sapiens/${VEP_VERSION}_GRCh38/.mva_complete"
chk "ClinVar（02 clinvar）" test -s "$CLINVAR"
chk "Exomiser ${EXOMISER_VERSION}（02 exomiser）" test -s "$EXOMISER_DIR/exomiser-cli-$EXOMISER_VERSION.jar"
for p in hg38 phenotype; do
  chk "Exomiser 資料 ${EXOMISER_DATA_VERSION}_${p}（02 exomiser）" test -f "$EXOMISER_DATA/${EXOMISER_DATA_VERSION}_$p/.mva_complete"
done
chk "MVA 基因組 BED（02 panel）" test -s "$PANEL_BED"
chk "HPO 本體 hp.obo（02 hpo）" test -s "$HPO_OBO"
chk "官方評分程式 evaluation.py（02 scorer）" test -s "$SCORER"

say "進度（${WORK}）"
for f in qc/stats.txt sex.txt norm.vcf.gz panel.vep.tsv panel.spliceai.tsv "exomiser/${PROBAND}_exome.genes.tsv" \
         "exomiser/${PROBAND}_introns.genes.tsv" trackB_pairs.tsv; do
  chk "$f" test -s "$WORK/$f"
done
if rules_done; then
  ok "排序規則已 commit（事先登錄）：$(git -C "$CONF" log -1 --format='%h %ci' -- ranking_rules.yml)"
else
  warn "排序規則還沒 commit：跑 06、07 之前要先 commit（見 README）"
fi
if ls "$RESULTS"/*.csv >/dev/null 2>&1; then ok "提交草稿：$(ls "$RESULTS"/*.csv | head -3 | tr '\n' ' ')"; else warn "還沒有提交草稿 CSV（08）"; fi
echo
