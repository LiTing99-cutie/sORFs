##### 0、准备文件 #####
# 20250910_c8_protein_map已经运行过
# bash run.prepare.20250912.sh

#####  2、所有样本 ##### 
# 更改参数，重新运行【先assign肽段，再进行过滤】【以这个为准】
# 20251010，添加source信息，重新运行
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
bash run_as.sh &> ../../log/run_as.log
## 追溯SNP
bash run_snp.sh &> ../../log/run_snp.log
## 查看SNP和质谱肽段的重叠
python3 snp_cov.py \
    ../../processed/variant_trace/snp_annotated.tsv \
    ../../results/all_samples.assignments.unique.post.tsv \
    ../../results/orfs_merged_final.tsv \
    ../../processed/variant_trace/manual_results.tsv
## report上一条
python3 gene_report.py ../../processed/variant_trace/manual_results.tsv
