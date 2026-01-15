# 输入文件路径
ORF_GTF="$1"
OUTPUT_PREFIX="$2"
VCF_FILE="/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/hom.vcf.gz"
# 参考基因组文件（用于突变类型判断，需要您提供路径）
REF_GENOME="/home/user/data/lit/database/public/genome/hg38/hg38.fa"  # 请修改为实际路径
python orf_snp.v2.20251121.py \
  --write-overlap $OUTPUT_PREFIX/snp_overlap.tsv \
  --orf-gtf $ORF_GTF \
  --ref-genome $REF_GENOME \
  --output $OUTPUT_PREFIX/snp_annotated.tsv \
  --vcf $VCF_FILE \
  --sample human_brain_21pcw