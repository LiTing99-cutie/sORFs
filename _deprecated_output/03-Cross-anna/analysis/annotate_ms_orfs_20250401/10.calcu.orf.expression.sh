mkdir -p output/S10
source /home/user/data2/lit/bin/lit_utils.sh 
define_annotation_gencode_v41_human
script=/home/user/data3/lit/project/sORFs/03-Cross-anna/Uni.gen.genepred.i_sorf_id.py
translate_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/S3.0c.Uni.translate_gtf.v2.20250325.sh
id_path=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S7/combined_id.txt
merged_bam_path=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S5/merged_bam
# 得到gtf
cd output/S10
python $script $id_path $gpe_15 gpe
genePredToGtf file gpe gtf
bash $translate_script gtf $fa
# 根据链特异性来计算表达量
## featureCounts的链特异性参数-s需要根据特定的文库来修改
gtf=$PWD/gtf
featureCounts -s 1 -T 10 -a $gtf -o rna_counts_1.txt $merged_bam_path/merged.1.bam
featureCounts -s 2 -T 10 -a $gtf -o rna_counts_2.txt $merged_bam_path/merged.2.bam
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o rna_counts_3.txt $merged_bam_path/merged.3.bam
featureCounts -T 10 -a $gtf -o rna_counts_4.txt $merged_bam_path/merged.4.bam