# 20241212修改脚本让一些参数可变
# 20241218修改卡的阈值
# 将每个软件输出的结果统一
# Usage: bash S3.0.Uni.Merge_res_v1.sh mouse_brain_output_20241011/*/*/output/Ribo*ORFs*
set -eo pipefail

work_path=$1
# 例如0.05或者0.1
pvalue_cutoff=$2
# 0或者1
only_cano_Scodon=$3
# custom或者youden或者fixed
riborf_cutoff=$4

cd $work_path
# 把之前的PRICE的结果软链接到这里，根据输出文件夹中的情况调整
# 20250128 注释这一行
# [ -d PRICE ] || ln -s ../Ribo-ORFs/PRICE/ ./
mkdir -p merge/RibORF merge/PRICE merge/RiboCode
cd merge
sample=$(echo $PWD | awk -F'/' '{print $(NF-3)}')
echo $sample

project_path=/home/user/data3/lit/project/sORFs/01-ribo-seq
# 20250125修改genePred
genePred=/home/user/data3/lit/project/sORFs/01-ribo-seq/output/assembled_trans/stringtie_output/gencode.vM29.add_assemble.genePred.txt
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
translate_script=$project_path/S3.0c.Uni.translate_gtf.sh
price_script=$project_path/S3.0b.Uni.Generate_genepred_PRICE_v1.py
compare_script_path=$project_path/ref/Ribo-seq-Tool-Comparison-Scripts-v2.0/Scripts_for_RiboCode_Analysis
riborf_script=$project_path/S3.0a.Uni.Format_RibORF_v1.R
ribocode_script=$project_path/S3.0f.Uni.Filter_RiboCode.R

# 使用base中的python
source activate base

# PRICE
cd PRICE
output_file=../../PRICE/*.orfs.tsv
## 根据PRICE的输出和已有注释的gpe构建一个sORF的gpe
python $price_script $output_file $genePred $pvalue_cutoff $only_cano_Scodon nonCano.formatted.gpe
## 根据gpe得到fasta文件
genePredToGtf file nonCano.formatted.gpe nonCano.formatted.gtf
bash $translate_script nonCano.formatted.gtf $fa
mv prot.fa nonCano.fa
## 得到长度后根据长度过滤
seqkit seq -g -m 6 -M 150 nonCano.fa > nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab

# RiboCode
cd ../RiboCode
output_1=../../RiboCode/$sample.txt
## 过滤结果
Rscript $ribocode_script $output_1 $pvalue_cutoff $only_cano_Scodon nonCano.sorf.txt
## 得到fasta文件
python $compare_script_path/Generate_Fasta_RiboCode.py nonCano.sorf.txt nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab

# RibORF
cd ../RibORF
output_1=../../RibORF/repre.valid.pred.pvalue.parameters.txt
output_2=../../RibORF/repre.valid.ORF.genepred.txt
stat_cutoff=../../RibORF/stat.cutoff.txt
Rscript $riborf_script $output_1 $output_2 $stat_cutoff $riborf_cutoff $only_cano_Scodon nonCano.sorf.formatted.gpe
genePredToGtf file nonCano.sorf.formatted.gpe nonCano.sorf.formatted.gtf
bash $translate_script nonCano.sorf.formatted.gtf $fa
mv prot.fa nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
