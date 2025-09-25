## 2025-07-10
### 运行fastqc
fastqc.sh
### Step 1-5: 变异检测主流程
nohup bash human_brain_21pcw_preprocessing.sh &> ../log/human_brain_21pcw_preprocessing.20250711.log &
- Conda环境：biotools
- 输入文件：
  - ../rawdata/L1EJF1602305-p21_Gen2Seq.R1.raw.fastq.gz
  - ../rawdata/L1EJF1602305-p21_Gen2Seq.R2.raw.fastq.gz
  - 参考基因组：/home/user/data/lit/database/public/genome/hg38/hg38.fa
- 线程数：30
- 输出文件：

## 2025-07-12 前处理
nohup bash human_brain_21pcw_preprocessing_v2_20250712.sh &> ../log/human_brain_21pcw_preprocessing_v2_20250712.20250712.log &

## 2025-07-14 碱基质量分数重校准；鉴定变异
nohup bash human_brain_21pcw_variant_calling_v1_20250714.sh &> ../log/human_brain_21pcw_variant_calling_v1_20250714.20250714.log &
- 备注：R语言没有安装相关的绘图包，导致中途运行暂停，应该index所有生成的bam
  
## 2025-07-15 碱基质量分数重校准；鉴定变异
nohup bash human_brain_21pcw_variant_calling_v2_20250715.sh &> ../log/human_brain_21pcw_variant_calling_v2_20250715.20250715.log &
- 安装了R语言相关包，且index了recal.bam，继续运行
- 安装命令为
  - conda activate biotools
  - Rscript -e "install.packages('gsalib', repos='http://cran.r-project.org')"
  - Rscript -e "install.packages(c('gplots', 'ggplot2'), repos='http://cran.r-project.org')"
  
## 2025-07-17 变异过滤
nohup bash filter_vcf_v1_20250717.sh &> ../log/filter_vcf_v1_20250717.20250717.log &
nohup bash calculate.filtered.snps_indels.number.20250718.sh &> ../log/calculate.filtered.snps_indels.number.20250718.20250718.log &

## 2025-07-18 变异注释
nohup bash annotate_vcf_vep_v1_20250718.sh &> ../log/annotate_vcf_vep_v1_20250718.20250718.log &

## 2025-07-18 变异注释第二次
path=/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter/
gunzip -c $path/human_brain_21pcw_filtered_pass.vcf.gz > $path/human_brain_21pcw_filtered_pass.vcf
bcftools view -H $path/human_brain_21pcw_filtered_pass.vcf|head -n 1000 >  $path/subset_1000.vcf
nohup bash annotate_vcf_vep_v1_20250718.1.sh &> ../log/annotate_vcf_vep_v1_20250718.20250718.1.log &
- 备注：上个命令退出状态不为0，重新调整了参数