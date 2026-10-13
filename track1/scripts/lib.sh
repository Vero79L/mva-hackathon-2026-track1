# shellcheck shell=bash
# Track 1 共用函式。每支腳本第一行：. "$(dirname "$0")/lib.sh"
# 相容 macOS 內建的 bash 3.2 與 BSD 指令（不用 GNU 專屬參數），Linux 也能跑。

set -euo pipefail

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
T1="$(cd "$SCRIPTS/.." && pwd -P)"
REPO="$(cd "$T1/.." && pwd -P)"
CONF="$T1/config"

# 先讀 local.env（個人覆蓋值，不進 git），再讀 paths.env（預設值）
# shellcheck disable=SC1091
if [ -f "$CONF/local.env" ]; then . "$CONF/local.env"; fi
# shellcheck source=../config/paths.env
. "$CONF/paths.env"
export SCRIPTS T1 REPO CONF MVA_ROOT DATA WORK RES RESULTS REF REF_PLAIN VEP_VERSION VEP_CACHE CLINVAR GTF \
  PANEL_FLANK PANEL_BED HPO_OBO SCORER EXOMISER_VERSION EXOMISER_DATA_VERSION EXOMISER_DIR EXOMISER_DATA \
  HF_DATASET HF_SPACE PROBAND THREADS VEP_FORK JAVA_XMX SPLICEAI_DIST

MAIN_CONTIG_RE='^(chr)?([0-9]+|X|Y|M|MT)$'

# ---------- 訊息 ----------
say()  { printf '\n==> %s\n' "$*"; }
ok()   { printf '  ✓ %s\n' "$*"; }
info() { printf '    %s\n' "$*"; }
warn() { printf '  ! %s\n' "$*" >&2; }
die()  { printf '\n✗ %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "找不到指令 $1。${2:-先跑 01_setup_mac.sh}"; }
is_mac() { [ "$(uname -s)" = Darwin ]; }

# ---------- 執行紀錄：畫面同時寫進 $WORK/logs（log 可能含變異資料，所以放 work，不放 repo） ----------
start_log() {
  mkdir -p "$WORK/logs"
  LOGFILE="$WORK/logs/$1_$(date +%Y%m%d-%H%M%S).log"
  exec > >(tee -a "$LOGFILE") 2>&1
  say "$1 開始於 $(date '+%Y-%m-%d %H:%M:%S')（紀錄：${LOGFILE}）"
}

# 長時間步驟：用 caffeinate 重新執行自己，避免 Mac 中途睡著（只在 macOS 有效）
stay_awake() {
  if [ -z "${MVA_CAFFEINATED:-}" ] && command -v caffeinate >/dev/null 2>&1; then
    export MVA_CAFFEINATED=1
    exec caffeinate -i "$BASH" "$0" "$@"
  fi
  export MVA_CAFFEINATED=1
}

# ---------- conda ----------
conda_base() {
  local d
  if [ -n "${CONDA_BASE:-}" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then echo "$CONDA_BASE"; return 0; fi
  for d in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
           /opt/homebrew/Caskroom/miniforge/base /usr/local/Caskroom/miniforge/base /opt/miniconda3 /opt/anaconda3; do
    if [ -f "$d/etc/profile.d/conda.sh" ]; then echo "$d"; return 0; fi
  done
  if command -v conda >/dev/null 2>&1; then conda info --base 2>/dev/null && return 0; fi
  return 1
}

has_env() { local b; b="$(conda_base)" || return 1; [ -d "$b/envs/$1" ]; }

conda_on() {   # conda_on <環境名稱>
  local b
  b="$(conda_base)" || die "找不到 conda：先跑 01_setup_mac.sh"
  set +u       # conda 的啟動腳本會讀未設定的變數
  # shellcheck disable=SC1091
  . "$b/etc/profile.d/conda.sh"
  if ! conda activate "$1"; then set -u; die "conda 環境 $1 不存在或無法啟動：先跑 01_setup_mac.sh"; fi
  set -u
}

# ---------- 檔案與磁碟 ----------
# fresh <輸出> <輸入...>：輸出存在、非空、而且比所有輸入新 → 回傳 0（可以跳過這一步）。FORCE=1 一律重做。
fresh() {
  local out="$1" f
  shift
  [ "${FORCE:-0}" = 1 ] && return 1
  [ -s "$out" ] || return 1
  for f in "$@"; do
    if [ -e "$f" ] && [ "$f" -nt "$out" ]; then return 1; fi
  done
  return 0
}

free_gb() {   # free_gb <路徑>：所在磁碟的可用空間（GB）；路徑還不存在時看最近一層存在的上層資料夾
  local d="$1"
  while [ ! -d "$d" ]; do d="$(dirname "$d")"; done
  df -Pk "$d" 2>/dev/null | awk 'NR==2 {printf "%d", $4/1048576}'
}

need_space() {   # need_space <路徑> <GB>
  local d="$1" f
  f="$(free_gb "$d")"
  [ "${f:-0}" -ge "$2" ] || die "$d 所在的磁碟只剩 ${f:-?} GB，這一步需要約 $2 GB。先清出空間，或在 config/local.env 把 RES 指到外接 SSD。"
}

# fetch <網址> <輸出檔>：可續傳；下載完成、大小核對過才改成正式檔名，所以「檔案存在」就代表完整
fetch() {
  local url="$1" out="$2" want have
  # 伺服器回 200 而且有給 Content-Length 才拿來核對。拿不到（離線、HEAD 被擋、錯誤頁）就不核對，
  # 也不能讓腳本在這裡無聲結束（set -e + pipefail），所以 curl 失敗要吞掉。
  want="$( (curl -sIL --max-time 60 "$url" 2>/dev/null || true) | tr -d '\r' \
    | awk '/^HTTP\//{code=$2; n=""} tolower($1)=="content-length:"{n=$2} END{if (code == 200) print n}')"
  if [ -s "$out" ]; then
    have="$(wc -c < "$out" | tr -d ' ')"
    if [ -z "$want" ] || [ "$have" = "$want" ]; then ok "已存在：$(basename "$out")"; return 0; fi
    # 大小不同：可能上次沒抓完，也可能伺服器換了新版（hp.obo 是 latest 網址）。接續會把兩個版本
    # 接成一個壞檔，所以從頭重抓到 .part；抓完、大小對了才取代舊檔。
    warn "$(basename "$out") 有 $have bytes，伺服器上是 $want bytes：從頭重新下載（抓完才取代舊檔）"
    rm -f "$out.part"
  fi
  mkdir -p "$(dirname "$out")"
  info "下載 $(basename "$out")"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$out.part" "$url" \
    || die "下載失敗：${url}（重跑同一支腳本會從中斷處續傳）"
  if [ -n "$want" ]; then
    have="$(wc -c < "$out.part" | tr -d ' ')"
    [ "$have" = "$want" ] || die "下載不完整：$(basename "$out") $have / $want bytes（重跑同一支腳本會續傳）"
  fi
  mv "$out.part" "$out"
}

unzip_big() {   # unzip_big <zip> <目的資料夾>：macOS 用 ditto（超過 4 GB 的 zip 也能解）
  mkdir -p "$2"
  if command -v ditto >/dev/null 2>&1; then ditto -x -k "$1" "$2"; else unzip -q -o "$1" -d "$2"; fi
}

sha256() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1; else sha256sum "$1" | cut -d' ' -f1; fi; }

# ---------- VCF ----------
SV_NAME_RE='(^|[/._-])(sv|svs|cnv|cnvs|manta|delly|lumpy|gridss|cnvnator|canvas|gcnv|dragen\.cnv|str|expansionhunter|mito)([/._-]|$)'

# find_vcf：主 VCF 的路徑。順序：環境變數 VCF → $DATA 裡唯一的小變異 VCF → 檔名 WGS_*.vcf.gz
find_vcf() {
  local hits n w
  if [ -n "${VCF:-}" ]; then
    [ -f "$VCF" ] || die "VCF=$VCF 不存在（檢查 config/local.env）"
    echo "$VCF"; return 0
  fi
  [ -d "$DATA" ] || die "資料夾 $DATA 不存在：先跑 03_fetch_data.sh"
  hits="$(find "$DATA" -type f -name '*.vcf.gz' 2>/dev/null | grep -viE "$SV_NAME_RE" | grep -vE '\.g\.vcf\.gz$|\.genome\.vcf\.gz$' | sort || true)"
  n="$(printf '%s\n' "$hits" | grep -c . || true)"
  if [ "$n" = 1 ]; then echo "$hits"; return 0; fi
  w="$(printf '%s\n' "$hits" | grep -E '/WGS_[^/]*\.vcf\.gz$' || true)"
  if [ "$(printf '%s\n' "$w" | grep -c . || true)" = 1 ]; then echo "$w"; return 0; fi
  if [ "$n" = 0 ]; then die "在 $DATA 找不到主 VCF：先跑 03_fetch_data.sh，或在 config/local.env 設定 VCF=完整路徑"; fi
  die "在 $DATA 找到不只一個 VCF，請在 config/local.env 指定 VCF=完整路徑：
$hits"
}

# sv_vcfs：資料夾裡看起來是 SV / CNV 的 VCF（每行一個）
sv_vcfs() { find "$DATA" -type f -name '*.vcf.gz' 2>/dev/null | grep -iE "$SV_NAME_RE" | sort || true; }

ensure_index() { [ -f "$1.tbi" ] || [ -f "$1.csi" ] || bcftools index -t "$1"; }

# main_contigs <vcf>：VCF 檔頭和參考序列都有的主要染色體（1–22、X、Y、M），逗號分隔
main_contigs() {
  bcftools view -h "$1" | sed -nE 's/^##contig=<ID=([^,>]+).*/\1/p' | grep -E "$MAIN_CONTIG_RE" > "$WORK/.contigs_vcf" || true
  cut -f1 "$REF.fai" | grep -Fx -f "$WORK/.contigs_vcf" | paste -sd, - || true
}

# ---------- 事先登錄檢查：ranking_rules.yml 必須已 commit、而且沒有未 commit 的修改 ----------
check_preregistered() {
  local h
  need git
  git -C "$CONF" rev-parse --git-dir >/dev/null 2>&1 || die "$REPO 不是 git repo，無法證明規則是事先登錄的。"
  h="$(git -C "$CONF" log -1 --format='%h %ci' -- ranking_rules.yml 2>/dev/null || true)"
  [ -n "$h" ] || die "排序規則還沒 commit（事先登錄）。先在終端機執行：
    cd \"$REPO\" && git add track1 && git commit -m \"Track 1: pipeline + pre-registered ranking rules v1\"
  然後再跑一次。"
  git -C "$CONF" diff --quiet HEAD -- ranking_rules.yml \
    || die "ranking_rules.yml 有還沒 commit 的修改（git diff 看改了什麼）。規則要先 commit 才能跑排序；改規則請把 version +1。"
  RULES_COMMIT="$h"
  export RULES_COMMIT
  ok "排序規則已事先登錄：commit $h"
}

rules_version() { sed -nE 's/^version:[[:space:]]*([0-9]+).*/\1/p' "$CONF/ranking_rules.yml" | head -1; }

# 執行紀錄（給方法報告用）：時間、步驟、規則 commit / 版本 / 雜湊
manifest() {   # manifest <步驟> <說明>
  mkdir -p "$RESULTS"
  local f="$RESULTS/run_manifest.tsv"
  [ -s "$f" ] || printf 'time\tstep\trules_commit\trules_version\trules_sha256\tnote\n' > "$f"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1" "${RULES_COMMIT:-}" \
    "$(rules_version)" "$(sha256 "$CONF/ranking_rules.yml")" "${2:-}" >> "$f"
}
