#!/usr/bin/env bash
# 用合成資料測試 Track 1 的程式邏輯（不需要比賽資料、VEP、Exomiser、網路，約 1 分鐘；共 45 個測試）。
# 非整倍體掃描的 3 個測試需要 bcftools + bgzip：mva 環境在就自動用它的，找不到就跳過（skip）。
# 用法：bash track1/tests/run_tests.sh
# 需要 pandas + PyYAML：有 conda 的 mva 環境就自動用它。
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd -P)"
py=python3
if ! "$py" -c 'import pandas, yaml' >/dev/null 2>&1; then
  for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
    if [ -x "$b/envs/mva/bin/python" ]; then py="$b/envs/mva/bin/python"; break; fi
  done
fi
"$py" -c 'import pandas, yaml' >/dev/null 2>&1 || { echo "需要 pandas 和 PyYAML：先跑 01_setup_mac.sh（或 conda activate mva）" >&2; exit 1; }
for b in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
  if [ -x "$b/envs/mva/bin/bcftools" ]; then PATH="${PATH}:$b/envs/mva/bin"; export PATH; break; fi
done
cd "$here"
exec "$py" -m unittest -v test_track1 test_phase2 test_benchmark
