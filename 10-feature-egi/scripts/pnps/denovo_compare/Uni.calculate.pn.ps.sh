##### 第四步、计算同义突变和非同义突变的数量 ######
vcf=$1
out_dir=$2
CDS_bed=$3
script=/home/user/data3/lit/project/sORFs/10-feature-egi/scripts/pnps/count_orf_syn_nonsyn_strict.py
# ====== 中间和输出路径 ======
tmp_dir=${out_dir}/tmp_strict
mkdir -p "$tmp_dir"

variants_bed=${tmp_dir}/variants.info.bed
intersect_tsv=${tmp_dir}/variants_in_orfs.strict.tsv
out_tsv=${out_dir}/orf_syn_nonsyn.strict.tsv

echo "[1/2] 从 VCF 提取变异位置 + REF/ALT/INFO ..."
# A: chr  start(0-based)  end(1-based)  REF  ALT  INFO
grep -v "^#" "$vcf" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2,$4,$5,$8}' \
  > "$variants_bed"

echo "[2/2] 与 ORF CDS bed 取交集 ..."
# B: chr start end ORF_ID ...
# 输出：1-6 是变异信息，7- 为 ORF bed 的列（第10列是 ORF_ID）
bedtools intersect -a "$variants_bed" -b "$CDS_bed" -wa -wb \
  > "$intersect_tsv"

echo "交集文件：$intersect_tsv"
echo "接下来运行 Python 严格统计脚本："
echo "  python count_orf_syn_nonsyn_strict.py $intersect_tsv $CDS_bed $out_tsv"
python $script $intersect_tsv $CDS_bed $out_tsv