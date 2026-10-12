#!/usr/bin/env bash
# 01 安裝工具鏈（macOS）：Xcode 命令列工具、Homebrew、Java 21、Miniforge、三個 conda 環境（mva / vep / spliceai）、
# 版本紀錄、Time Machine 排除、git 防護。可以重複執行：已經裝好的會跳過。約 30–60 分鐘。
# 用法：bash track1/scripts/01_setup_mac.sh [--no-spliceai]
. "$(dirname "$0")/lib.sh"

with_spliceai=1
for arg in "$@"; do
  case "$arg" in
    --no-spliceai) with_spliceai=0 ;;
    *) die "不認識的參數：${arg}（只有 --no-spliceai）" ;;
  esac
done
is_mac || die "這支腳本是 macOS 版。Linux 請用 conda 裝同樣版本的 bcftools、ensembl-vep=${VEP_VERSION}、Java 21。"
start_log 01_setup

say "1/8 Xcode 命令列工具"
if xcode-select -p >/dev/null 2>&1; then ok "已安裝"; else
  xcode-select --install || true
  die "跳出的視窗按「安裝」，裝完再跑一次這支腳本"
fi

say "2/8 Homebrew"
BREW=""
for b in /opt/homebrew/bin/brew /usr/local/bin/brew; do
  if [ -x "$b" ]; then BREW="$b"; break; fi
done
if [ -z "$BREW" ]; then
  die "還沒有 Homebrew。把下面這行貼進終端機（會要求輸入 Mac 密碼），裝完再跑一次這支腳本：
  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
fi
eval "$("$BREW" shellenv)"
grep -q 'brew shellenv' "$HOME/.zprofile" 2>/dev/null || printf '\neval "$(%s shellenv)"\n' "$BREW" >> "$HOME/.zprofile"
ok "$(brew --version 2>/dev/null | head -1 || true)"

say "3/8 Homebrew 套件與 Java 21"
brew list --formula gh >/dev/null 2>&1 || brew install gh
if /usr/libexec/java_home -v 21+ >/dev/null 2>&1; then
  ok "Java：$(java -version 2>&1 | head -1 || true)"
else
  info "安裝 Temurin 21（會要求輸入 Mac 密碼）"
  brew install --cask temurin@21
fi

say "4/8 Miniforge（conda / mamba）"
if base="$(conda_base)"; then
  ok "已有 conda：$base"
else
  inst="$HOME/Miniforge3-installer.sh"
  curl -fL -o "$inst" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"
  bash "$inst" -b -p "$HOME/miniforge3"
  rm -f "$inst"
  base="$HOME/miniforge3"
  "$base/bin/conda" init zsh >/dev/null && ok "已安裝並加進 ~/.zshrc（新開的終端機就能用 conda activate）"
fi
SOLVER="$base/bin/mamba"
[ -x "$SOLVER" ] || SOLVER="$base/bin/conda"

ensure_rosetta() {
  [ "$(uname -m)" = arm64 ] || return 0
  if /usr/bin/pgrep -q oahd 2>/dev/null; then return 0; fi
  info "安裝 Rosetta 2（Intel 版套件需要）"
  softwareupdate --install-rosetta --agree-to-license \
    || die "Rosetta 安裝失敗：在終端機執行  sudo softwareupdate --install-rosetta --agree-to-license  之後再跑一次"
}

# create_env <名稱> <套件...>：先試 Apple Silicon 原生版；沒有原生版本時改建 Intel 版，由 Rosetta 執行
create_env() {
  local name="$1"
  shift
  if has_env "$name"; then ok "conda 環境 $name 已存在"; return 0; fi
  info "建立 conda 環境 ${name}：$*"
  if "$SOLVER" create -y -n "$name" -c conda-forge -c bioconda --override-channels --strict-channel-priority "$@"; then
    return 0
  fi
  [ "$(uname -m)" = arm64 ] || return 1
  warn "$name 沒有 Apple Silicon 原生版本，改建 Intel 版（Rosetta 執行）"
  ensure_rosetta
  rm -rf "${base:?}/envs/$name"
  CONDA_SUBDIR=osx-64 "$SOLVER" create -y -n "$name" -c conda-forge -c bioconda --override-channels \
    --strict-channel-priority "$@" || return 1
  "$base/bin/conda" config --file "$base/envs/$name/.condarc" --set subdir osx-64   # 之後在這個環境裝套件也用 Intel 版
}

say "5/8 conda 環境"
create_env mva python=3.11 "bcftools>=1.20" "htslib>=1.20" "samtools>=1.20" pandas pyyaml openpyxl pypdf requests \
  || die "mva 環境建立失敗：看上面的錯誤訊息"
"$base/envs/mva/bin/python" -m pip install -q -U "huggingface_hub[hf_xet]"
ok "mva：$("$base/envs/mva/bin/bcftools" --version 2>/dev/null | head -1 || true)"

create_env vep "ensembl-vep=$VEP_VERSION" \
  || die "vep 環境建立失敗。用  CONDA_SUBDIR=osx-64 mamba search -c bioconda ensembl-vep  查有哪些版本；
  若沒有 ${VEP_VERSION}，在 config/local.env 設 VEP_VERSION=<最新版號>（cache 會跟著下載同一版）再跑一次。"
if PATH="$base/envs/vep/bin:$PATH" vep --help 2>/dev/null | grep -E "ensembl-vep +: +$VEP_VERSION" >/dev/null; then
  ok "vep：ensembl-vep $VEP_VERSION"
else
  warn "vep --help 沒看到版本 ${VEP_VERSION}，請手動確認：conda activate vep && vep --help | head"
fi

if [ "$with_spliceai" = 1 ]; then
  if create_env spliceai python=3.10 pip \
     && "$base/envs/spliceai/bin/python" -m pip install -q "tensorflow==2.15.*" "spliceai==1.3.1" \
     && "$base/envs/spliceai/bin/spliceai" -h >/dev/null 2>&1; then
    ok "spliceai：已安裝（TensorFlow 2.15）"
  else
    warn "SpliceAI 沒裝起來。第一階段可以先跳過：排序會把 SpliceAI 分數當 0，第二階段再補。"
  fi
else
  info "略過 SpliceAI（--no-spliceai）"
fi

say "6/8 版本紀錄（$T1/env，進 repo 給報告用）"
mkdir -p "$T1/env"
for e in mva vep spliceai; do
  if has_env "$e"; then "$base/bin/conda" env export -n "$e" --no-builds 2>/dev/null | grep -v '^prefix:' > "$T1/env/$e.yml" || true; fi
done
{ sw_vers; uname -m; } > "$T1/env/macos.txt"
java -version 2>&1 | head -1 > "$T1/env/java.txt" || true
{
  "$base/envs/mva/bin/bcftools" --version | head -1
  "$base/envs/mva/bin/samtools" --version | head -1
  "$base/envs/mva/bin/python" --version
  "$base/envs/mva/bin/python" -c 'import huggingface_hub as h, pandas as p; print("huggingface_hub", h.__version__, "| pandas", p.__version__)'
  PATH="$base/envs/vep/bin:$PATH" vep --help 2>/dev/null | grep -E '^ +ensembl' | head -5
} > "$T1/env/tools.txt" 2>/dev/null || true
ok "$(ls "$T1/env" | tr '\n' ' ')"

say "7/8 資料夾與 Time Machine 排除"
mkdir -p "$DATA" "$WORK" "$RES" "$RESULTS"
for d in "$DATA" "$MVA_ROOT/work" "$WORK"; do
  [ -d "$d" ] || continue
  if tmutil isexcluded "$d" 2>/dev/null | grep Excluded >/dev/null; then ok "Time Machine 已排除 $d"
  elif tmutil addexclusion "$d"; then ok "Time Machine 排除 $d"
  else warn "無法排除 ${d}：到「系統設定 → 一般 → Time Machine → 選項」手動加入"; fi
done
info "Spotlight 只能手動：系統設定 → Spotlight →「搜尋隱私權」，把 $WORK 加進去"

say "8/8 git 防護"
if git -C "$REPO" rev-parse --git-dir >/dev/null 2>&1; then
  hook="$(git -C "$REPO" rev-parse --git-path hooks/pre-commit)"
  case "$hook" in /*) ;; *) hook="$REPO/$hook" ;; esac
  mkdir -p "$(dirname "$hook")"
  if [ ! -e "$hook" ]; then
    cp "$SCRIPTS/git-pre-commit.sh" "$hook" && chmod +x "$hook"
    ok "已安裝 pre-commit 防護：擋基因體檔案、超過 5 MB 的檔案、gated-data / work 底下的檔案"
  elif grep -q 'mva-track1-guard' "$hook"; then
    ok "pre-commit 防護已安裝"
  else
    warn "已經有別的 pre-commit hook，沒有覆蓋。要加防護：把 track1/scripts/git-pre-commit.sh 的內容併進 $hook"
  fi
else
  warn "$REPO 不是 git repo：git init 之後再跑一次這支腳本"
fi

conda_on mva
if hf auth whoami >/dev/null 2>&1; then
  ok "Hugging Face 已登入"
else
  warn "Hugging Face 還沒登入：執行  conda activate mva && hf auth login  （貼 Read token；問要不要存成 git credential 選 n）"
fi
say "完成。下一步：bash track1/scripts/02_fetch_resources.sh（大檔案多，建議晚上接電源掛著下載）"
