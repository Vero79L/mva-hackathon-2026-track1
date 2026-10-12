#!/usr/bin/env bash
# mva-track1-guard：git pre-commit 防護（01_setup_mac.sh 會複製到 .git/hooks/pre-commit）
# 擋下：基因體檔案、gated-data / work 底下的檔案、超過 5 MB 的檔案、看起來是 Exomiser / 路線 B 結果的表格。
# 真的需要提交大檔時（例如圖片），用 git commit --no-verify 跳過，並自己確認內容。
blocked=0
files="$(git diff --cached --name-only --diff-filter=ACM)"
[ -z "$files" ] && exit 0

hits="$(printf '%s\n' "$files" | grep -Ei '\.(vcf|bcf|bam|cram|fastq|fq|gvcf)(\.gz)?$|\.(tbi|csi|bai|crai)$|(^|/)(gated-data|work)/|\.(genes|variants)\.tsv$|trackB_.*\.tsv$|panel\.(vep|spliceai)\.tsv$|phenopacket\.ya?ml$|hpo_terms\.tsv$' || true)"
if [ -n "$hits" ]; then
  echo "BLOCKED（mva-track1-guard）：這些檔案可能含病人資料，不能 commit：" >&2
  printf '%s\n' "$hits" | sed 's/^/  /' >&2
  blocked=1
fi

printf '%s\n' "$files" | while IFS= read -r f; do
  if [ -f "$f" ] && [ "$(wc -c < "$f" | tr -d ' ')" -gt 5000000 ]; then
    echo "BLOCKED（mva-track1-guard）：$f 超過 5 MB" >&2
    exit 1
  fi
done || blocked=1

exit "$blocked"
