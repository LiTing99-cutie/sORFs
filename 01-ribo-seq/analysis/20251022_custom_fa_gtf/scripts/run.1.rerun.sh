# 0、整理clean reads的文件list
proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis
output_path=$proj_path/20251022_custom_fa_gtf/processed
cat <(find $proj_path/20250722_formal_data_run/ -name "*trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz") \
 <(find $proj_path/20250813_demo_data_analysis/ -name "*trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz") > ../processed/clean_fastq.lst

# 1、从clean reads开始比对【更改了参考基因组】
# echo -e "***Mapping at $(date '+%Y-%m-%d %H:%M:%S')"
# bash mapping.20251022.sh &> ../log/mapping.log

# 2、提取特定长度的reads【后续读端长度分布、基因组位置分布以及使用ribotish计算offset都是基于此进行计算】
# echo -e "***Extract reads at $(date '+%Y-%m-%d %H:%M:%S')"
# bash run.extract.25-34.length.reads.20250725.sh $output_path/alignment $output_path/filtered_bam

# 3、bam合并
echo -e "***Merge bam at $(date '+%Y-%m-%d %H:%M:%S')"
mergd_bam_out_path=$(realpath ../processed/merged_bam)
mkdir -p $mergd_bam_out_path
## 合并所有的PRF bam【25-34长度】
# samtools merge -@ 30 -f -o $mergd_bam_out_path/merged.bam \
#   $output_path/filtered_bam/*.bam && samtools index -@ 30 $mergd_bam_out_path/merged.bam
## 合并所有的toTranscriptome.out.bam以及sortedByCoord.out.bam供ORF鉴定分析
# out=$mergd_bam_out_path/Aligned.toTranscriptome.out.bam
# # 需要更改
# samtools merge -@ 30 -f -o $out \
#   $(find $output_path/alignment/ -name "*Aligned.toTranscriptome.out.bam")
# out=$mergd_bam_out_path/Aligned.sortedByCoord.out.bam
# samtools merge -@ 30 -f -o $out \
#   $(find $output_path/alignment/ -name "*Aligned.sortedByCoord.out.bam")

# 4、鉴定ORF
echo -e "***Call ORFs at $(date '+%Y-%m-%d %H:%M:%S')"
bash call-orfs/call.orfs.20251024.sh $(realpath ../processed/orfs) \
  $mergd_bam_out_path/Aligned.toTranscriptome.out.bam \
  $mergd_bam_out_path/Aligned.sortedByCoord.out.bam  &> ../log/call.orfs.log

# 5、整理鉴定出的ORF
echo -e "***Organize ORFs at $(date '+%Y-%m-%d %H:%M:%S')"
## 1：工作路径；2：PRICE软件的cutoff；3：是否只需要经典的起始密码子；4：riborf_cutoff如何选取
## 5：是否只需要经典的ORF 6：工作路径下的输出文件夹名
# 0或者1
bash organize_res/2.3.Uni.Organize_res_v3.20251016.sh ../processed/orfs \
  0.05 0 custom 0 organized
bash organize_res/3.1.Uni.Merge_Filter_Annotate.v2.20250424.sh ../processed/orfs/organized ../processed/orfs/merged