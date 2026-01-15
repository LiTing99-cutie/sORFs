#!/bin/bash
set -euo pipefail
source activate biotools
# ====== 输入文件 ======
vcf=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/pnps/orf_regions/hg38_custom.eff.vcf
orf_bed=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/intergenic.lnc.noncano.cano.orfs.bed

# ====== 中间和输出路径 ======
base_out=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/pnps
tmp_dir=${base_out}/tmp_strict
mkdir -p "$tmp_dir"

variants_bed=${tmp_dir}/variants.info.bed
intersect_tsv=${tmp_dir}/variants_in_orfs.strict.tsv
out_tsv=${base_out}/orf_syn_nonsyn.strict.tsv

echo "[1/2] 从 VCF 提取变异位置 + REF/ALT/INFO ..."
# A: chr  start(0-based)  end(1-based)  REF  ALT  INFO
bcftools view -H "$vcf" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2,$4,$5,$8}' \
  > "$variants_bed"

echo "[2/2] 与 ORF CDS bed 取交集 ..."
# B: chr start end ORF_ID ...
# 输出：1-6 是变异信息，7- 为 ORF bed 的列（第10列是 ORF_ID）
bedtools intersect -a "$variants_bed" -b "$orf_bed" -wa -wb \
  > "$intersect_tsv"

echo "交集文件：$intersect_tsv"
echo "接下来运行 Python 严格统计脚本："
echo "  python count_orf_syn_nonsyn_strict.py $intersect_tsv $orf_bed $out_tsv"
python pnps/count_orf_syn_nonsyn_strict.py $intersect_tsv $orf_bed $out_tsv