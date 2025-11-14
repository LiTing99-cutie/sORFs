# 输入文件路径
VCF_FILE="/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/hom.vcf.gz"
ORF_GTF="/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf"
ANNOTATION_GTF="/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf"
CLASSIFICATION_FILE="/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/classify/collapsed_classification.filtered_lite_classification.txt"
# 参考基因组文件（用于突变类型判断，需要您提供路径）
REF_GENOME="/home/user/data/lit/database/public/genome/hg38/hg38.fa"  # 请修改为实际路径
OUTPUT_PREFIX="/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/processed/variant_trace"
python orf_snp.py \
  --write-overlap $OUTPUT_PREFIX/snp_overlap.tsv \
  --orf-gtf $ORF_GTF \
  --ref-genome $REF_GENOME \
  --output $OUTPUT_PREFIX/snp_annotated.tsv \
  --vcf $VCF_FILE \
  --sample human_brain_21pcw