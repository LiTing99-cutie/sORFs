##### 0、准备文件 #####
# 20250910_c8_protein_map已经运行过
# bash run.prepare.20250912.sh

##### 1、合并public以及in house的质谱搜库结果 #####
path_1=../MS_res_from_Galaxy/public
path_2=../MS_res_from_Galaxy
path_3=../MS_res_from_Galaxy/merge_public_in_house
for filename in merged.peptide.tsv pep.orf.merged.txt peptide_intensity_IL.merged.tsv;do
head -n2 $path_1/$filename
head -n2 $path_2/$filename
cat $path_1/$filename <(tail -n +2 $path_2/$filename) > $path_3/$filename
done
#####  2、所有样本 ##### 
# 更改参数，重新运行【先assign肽段，再进行过滤】【以这个为准】
# 20251010，添加source信息，重新运行
# 20251119，修改了产生pep orf时追加而不是写入的问题，split_by_sample.v2.20251119.py
# 20251119，修改了新样本没有酶切后缀无法定义酶切规则的问题，count_theoretical_peptides.v1.20251119.py
## 如果有新的样本，建议使用这个脚本先小测一下再批量运行
## 打开看下run.all.sample.20250918.sh里面质谱相关的文件路径需要更改
## 酶切规则默认trypsin
# bash run.all.steps.20250918.v3.sh --sample 59GM1 --by-sample-dir ../processed/by_sample      --rpf-min 0 --ps-min 0 --min-c 0
nohup bash run.all.sample.20250918.sh \
  --by-sample-dir ../processed/by_sample/ \
  --logroot ../log/run_all_sample \
  --results-dir ../results/ \
  --nproc 8 \
  --rpf-min 0 --ps-min 0 --min-c 0 &>../log/run.all.sample.20251010.log &
## 添加起始位点和终止位点以及CDS序列
mkdir -p ../results/augment_orf_table
bash augment_orf_table.sh ../results/orfs_merged_final.tsv ../results/augment_orf_table
#####  3、查看kozak等pattern ##### 
#####  4、查看变异来源 ##### 
cd variant_trace
## 追溯可变剪切
bash run_as.sh ../../results/augment_orf_table/sub.gtf ../../processed/variant_trace/ &> ../../log/run_as.log
## 追溯SNP
bash run_snp.sh ../../results/augment_orf_table/sub.gtf ../../processed/variant_trace/ &> ../../log/run_snp.log
## 查看SNP和质谱肽段的重叠
python3 snp_cov.py \
    ../../processed/variant_trace/snp_annotated.tsv \
    ../../results/all_samples.assignments.unique.post.tsv \
    ../../results/orfs_merged_final.tsv \
    ../../processed/variant_trace/manual_results.tsv
## report上一条
python3 gene_report.py ../../processed/variant_trace/manual_results.tsv

#####  5、可视化 ##### 