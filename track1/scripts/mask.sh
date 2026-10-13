# shellcheck shell=bash
# mask.sh — STATUS.md 的遮罩引擎（CLAUDE.md 第 5 節）。遮罩一律由程式做，不靠 LLM。
# 獨立於 lib.sh：測試可以單獨 source。實際比對用 perl（macOS / Linux 都內建），
# 靜態規則寫死在這裡；執行時詞表（樣本名、gated-data 檔名）由 mva_mask_init 收集。

MVA_MASK_RULES="${MVA_MASK_RULES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)/config/ranking_rules.yml}"

# 基因名與別名：從 ranking_rules.yml 讀（tier_a、tier_b、aliases 的 key 與 value）；讀不到就用內建清單
mva_mask_genes() {
  local yml="$MVA_MASK_RULES" names=""
  if [ -f "$yml" ]; then
    names="$( { sed -nE 's/^ *tier_[ab]: *\[([^]]*)\].*/\1/p' "$yml" | tr ',' '\n'
                sed -nE 's/^ *aliases: *\{([^}]*)\}.*/\1/p' "$yml" | tr ',:' '\n'; } )"
    # shellcheck disable=SC2086
    names="$(printf '%s\n' $names | grep -E '^[A-Z][A-Z0-9]+$' | sort -u || true)"
  fi
  if [ -z "$names" ]; then
    names="BUB1 BUB1B BUB3 CASC5 CCDC84 CENATAC CEP192 CEP57 FAM178A KNL1 MAD1L1 MAD2L1BP SLF2 SMC5 TRIP13 ZWINT"
  fi
  # shellcheck disable=SC2086
  printf '%s\n' $names | paste -sd'|' -
}

# 執行時詞表：VCF header 的樣本名 + gated-data 底下所有檔名。只有 supervisor 會呼叫。
mva_mask_init() {
  local vcf words=""
  if type find_vcf >/dev/null 2>&1 && command -v bcftools >/dev/null 2>&1; then
    if vcf="$(find_vcf 2>/dev/null)"; then
      words="$(bcftools query -l "$vcf" 2>/dev/null || true)"
    fi
  fi
  if [ -n "${DATA:-}" ] && [ -d "$DATA" ]; then
    words="$words
$(find "$DATA" -type f -exec basename {} \; 2>/dev/null || true)"
  fi
  MVA_MASK_WORDS="$words"
  export MVA_MASK_WORDS
}

# mva_mask：stdin → stdout。順序：詞表 → gated-data 路徑 → HGVS → rsID → HPO → 座標
# → 等位基因 → 鹼基序列 → 基因型 → 基因名 → 性別 → 4 位以上數字（最後，前面規則已把特定格式收掉）。
mva_mask() {
  MVA_MASK_GENES="$(mva_mask_genes)" perl -pe '
    BEGIN {
      @mask_words = grep { length } split /\n/, ($ENV{MVA_MASK_WORDS} || "");
      $mask_genes = $ENV{MVA_MASK_GENES} || "";
    }
    for $w (@mask_words) { $q = quotemeta($w); s/$q/<masked>/g; }
    s{gated-data/\S+}{gated-data/<masked>}g;
    s/\b[cpgn]\.(?=[A-Za-z0-9*>+=_-]*[0-9])[A-Za-z0-9*>+=_.-]+/<hgvs>/g;
    s/\brs[0-9]+\b/<rsid>/g;
    s/\bHP:[0-9]+\b/<hpo>/g;
    s/\b(chr)?([0-9]{1,2}|X|Y|MT?):[0-9]+(-[0-9]+)?\b/<coord>/g;
    s/\b[ACGTN]+>[ACGTN]+\b/<allele>/g;
    s/\b[ACGTN]{4,}\b/<seq>/g;
    s{(?<![0-9.\/|])[0-9.][\/|][0-9.](?![0-9.\/|])}{<gt>}g;
    if (length $mask_genes) { s/\b($mask_genes)\b/<gene>/g; }
    s/\b(FE)?MALE\b/<sex>/gi;
    s/[0-9]{4,}/<num>/g;
  '
}
