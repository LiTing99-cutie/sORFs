##### 0、准备文件 #####
bash run.prepare.20250912.sh

##### 1、单个样本测试 #####
id_map=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt
bash id.convert.20250911.sh ../processed/pick_orf/pep.orf.txt \
                   "$id_map" \
                   ../processed/pick_orf/pep.orf.rePicked.txt
python pep_to_protein_assign.20250915.alt.col.py \
  --pep2prot ../processed/pick_orf/pep.orf.rePicked.txt \
  --isoform /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --orf_psites /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --nca /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --rpf_metric RPF_codon_coverage --rpf_min 0.2 \
  --ps_metric Psites_codon_coverage --ps_min 0.1 \
  --min_c 0.2 \
  --outdir ../processed/pep_assign

python pep_to_protein_assign.20250915.alt.col.py \
  --pep2prot ../processed/pick_orf/pep.orf.rePicked.txt \
  --isoform /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --orf_psites /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --nca /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --rpf_metric RPF_codon_coverage --rpf_min 0 \
  --ps_metric Psites_codon_coverage --ps_min 0 \
  --min_c 0 \
  --outdir ../processed/pep_assign

python make_orf_tables.py \
  --cand_filtered ../processed/pep_assign/candidates.filtered.tsv \
  --orf_seq_len   "$OUT/orf_seq_len.tsv" \
  --isoform_info  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --out_detail    ../processed/pep_assign/orf_detail.noCanonSharing.tsv \
  --out_orf       ../processed/pep_assign/orf_folded_counts.tsv

sample=21pcw_1_C8_T_T
python3 count_theoretical_peptides.py \
  --in ../processed/pep_assign/orf_folded_counts.tsv \
  --sample $sample \
  --out ../processed/pep_assign/orf_folded_counts_with_theo_pep.tsv

cat <(head -n1 ../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv) \
 <(grep 1_C8_T_T ../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv) > ../processed/quant/peptide_intensity_IL.tsv

python compute_ibaq_len_norm.py \
  --file1 ../processed/pep_assign/assignments.unique.post.tsv \
  --file2 ../processed/quant/peptide_intensity_IL.tsv \
  --file3 ../processed/pep_assign/orf_folded_counts.tsv \
  --file3-cols all \
  -o ../processed/quant/ibaq_b_with_total.tsv

python compute_ibaq_len_norm.v1.py \
  --file1 ../processed/pep_assign/assignments.unique.post.tsv \
  --file2 ../processed/quant/peptide_intensity_IL.tsv \
  --file3 ../processed/pep_assign/orf_folded_counts_with_theo_pep.tsv \
  --file3-cols all \
  -o ../processed/quant/ibaq_b_with_total.tsv

python3 iBAQ_A_B_correlation.py \
  --input ../processed/quant/ibaq_b_with_total.tsv \
  --log \
  --plot ../processed/quant/iBAQ_A_vs_B.png
# 以上已经整合到通用脚本run.all.20250912.sh中

# 合并定量的信息【所有样本运行完毕后再整合】
python merge_tables.py \
  --ibaq ../processed/quant/ibaq_b_with_total.tsv \
  --rpf  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --iso  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --rna  /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --gene_anno /home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt \
  --out  ../results/ibaq_orf_rpf_iso_rna.tsv

# python3 ../_deprecated_scripts/spearman_corr_single_sample.py \
#   --merged ../results/ibaq_orf_rpf_iso_rna.tsv \
#   --out ../results/spearman_results.tsv

python3 spearman_corr.v2.py \
  --merged ../results/ibaq_orf_rpf_iso_rna.tsv \
  --out ../results/spearman_results.tsv \
  --pairs A:RPF_RPKM A:iBAQ_B RPF_RPKM:iBAQ_B
#           group     var_x     var_y  spearman_rho     n note
# 0           all         A  RPF_RPKM      0.573695  1915     
# 1           all         A    iBAQ_B      0.302686  1543     
# 2           all  RPF_RPKM    iBAQ_B      0.273301  1543 
python3 spearman_corr.v2.py \
  --merged ../results/ibaq_orf_rpf_iso_rna.tsv \
  --out ../results/spearman_results.tsv \
  --pairs A:RPF_RPKM A:iBAQ_A RPF_RPKM:iBAQ_A
#           group     var_x     var_y  spearman_rho     n note
# 0           all         A  RPF_RPKM      0.573695  1915     
# 1           all         A    iBAQ_A      0.285343  1543     
# 2           all  RPF_RPKM    iBAQ_A      0.257315  1543 

# 单个样本看看时间【2-3min】
time bash run.all.20250912.sh 21pcw_1_C8_T_T

#####  2、所有样本 ##### 
# 更改参数，重新运行【先assign肽段，再进行过滤】【以这个为准】
# 20251010，添加source信息，重新运行
nohup bash run.all.sample.20250918.sh \
  --by-sample-dir ../processed/by_sample/1 \
  --logroot ../logs/1 \
  --results-dir ../results/1 \
  --nproc 8 \
  --rpf-min 0 --ps-min 0 --min-c 0 &>../log/1/run.all.sample.20251010.log &
## 添加起始位点和终止位点以及CDS序列
mkdir -p ../results/1/augment_orf_table
bash augment_orf_table.sh ../results/1/orfs_merged_final.tsv ../results/1/augment_orf_table
#####  3、sORF很多，但是kozak等pattern并不好 ##### 
# 以C8为例，查看肽段的来源
sample=21pcw_1_C8_T_T
# 31106
awk '$4=="'$sample'"' ../MS_res_from_Galaxy/merged.peptide.tsv|wc -l
# 26357
awk '$4=="'$sample'"' ../MS_res_from_Galaxy/merged.peptide.tsv|grep msfragger_closed| wc -l
# 9387
join -1 1 -2 1 <(sort -k1,1 ../processed/by_sample/1/21pcw_1_C8_T_T/pep_assign/assignments.unique.post.tsv) \
 <(awk '$4=="'$sample'"' ../MS_res_from_Galaxy/merged.peptide.tsv|sort -k1,1)|wc -l
# 7651
join -1 1 -2 1 <(sort -k1,1 ../processed/by_sample/1/21pcw_1_C8_T_T/pep_assign/assignments.unique.post.tsv) \
 <(awk '$4=="'$sample'"' ../MS_res_from_Galaxy/merged.peptide.tsv|sort -k1,1)|grep msfragger_closed|wc -l

awk -v OFS='\t' '{print ".",$1}' ../MS_res_from_Galaxy/ms_closed.unique.protein.txt > ../processed/tmp/ms_closed.unique.protein.txt
proj_path=/home/user/data3/lit/project/sORFs/
id_map=$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt
bash id.convert.20250911.sh ../processed/tmp/ms_closed.unique.protein.txt \
                   "$id_map" \
                   ../processed/tmp/pep.orf.rePicked.txt
# 17538
wc -l  ../processed/tmp/pep.orf.rePicked.txt
#####  4、在建模前，过滤掉那些和经典的ORF有重叠的ORF，顺便计算下in-frame的重叠的比例 #####
conda activate base
sep_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/sub.gtf
anno_gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
python3 orf_overlap_inframe.py --new-gtf $sep_gtf --ann-gtf $anno_gtf --out ../results/1/orf_overlap_inframe.txt
#####  5、计算periodicity ##### 
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250909_org_all_data/processed/all.offsetCorrected.merged.sorted.bam
gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/sub.gtf
mkdir -p ../results/1/peri/
python3 calcu_peri/psite_frame_stats.v1.py \
  --bam $bam \
  --annot $gtf \
  --format gtf \
  --key-attr gene_id \
  --out ../results/1/peri/psite_frame_stats.v1.tsv \
  --log-every 2000