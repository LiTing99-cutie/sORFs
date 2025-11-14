# 1、评估ribo-seq文库结果【不统计原始reads和去污染步骤中的reads】
proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis
output_path=$proj_path/20251022_custom_fa_gtf/processed
bash qual_assess/ribo.qual.assess.20251022.sh $output_path
bash qual_assess/ribo.qual.assess.20251024.raw.contam.reads.sh $output_path
mkdir -p ../figures/qual_assess
Rscript qual_assess/Uni.organize.v1.20250730.R $(realpath ../processed/qual_assess/) \
    $(realpath ../figures/qual_assess)
# 图片的美化
Rscript qual_assess/Uni.organize.v2.figure.sample.reorder.R $(realpath ../processed/qual_assess/) \
    $(realpath ../figures/qual_assess)
# 导出制作表格
{ head -n1 ../figures/qual_assess/stat.all.txt; \
tail -n +2 ../figures/qual_assess/stat.all.txt | sort -t '_' -k1,1 -k2,2n -k3,3n; } | \
cut -f 1-8,11 > ../figures/qual_assess/stat.all.org.txt

awk -F'\t' 'NR==1{print; next} {
  # 仅当第8列是数值时格式化；否则保持原样
  if ($8 ~ /^[-+]?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?$/) $8 = sprintf("%.2f", $8)
  print
}' OFS='\t' ../figures/qual_assess/stat.all.org.txt > ../figures/qual_assess/stat.all.org.1.txt

# 2、测序饱和分析
bash run.saturation.20251026.sh