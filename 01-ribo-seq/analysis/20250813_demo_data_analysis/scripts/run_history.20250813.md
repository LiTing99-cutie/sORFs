# 1一、qc，去除污染，比对
## 生成原始数据list文件
proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq
raw_data_path=$proj_path/rawdata/demo_20250813/org_rename
ls $raw_data_path/*.fq.gz > ../processed/raw_fastq_batch_1.lst
## qc去接头，去除污染，比对
processed_dir=$proj_path/analysis/20250813_demo_data_analysis/processed
nohup bash Uni.qc.mapping.using.star.v3.20250725.sh $processed_dir/raw_fastq_batch_1.lst $processed_dir/batch_1 &> ../log/Run.qc.mapping.using.star.v3.log &

# 二、提取特定长度的reads
nohup bash run.extract.25-34.length.reads.20250725.sh &
- 提取特异性比对的reads中符合特定长度的reads

# 四、评估ribo-seq文库结果
bash ribo.qual.assess.20250730.sh
mkdir -p ../figures/batch_1/qual_assess
Rscript Uni.organize.v1.20250730.R $(realpath ../results/batch_1/qual_assess) \
    $(realpath ../figures/batch_1/qual_assess)


# 三、测序饱和分析
## 合并所有的bam
mkdir -p ../processed/batch_1_2_merged
nohup samtools merge -@ 30 -f -o ../processed/batch_1_2_merged/batch_1_2_merged.bam \
  /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1/filtered_bam/*.bam \
  ../processed/batch_1/filtered_bam/*.bam  \
&& samtools index -@ 30 ../processed/batch_1_2_merged/batch_1_2_merged.bam &

## 饱和分析
### annotated gtf; allow overlap
bash run.saturation.20250815.sh
### annotated gtf; do not allow overlap
bash run.saturation.20250815.1.sh
### iso-seq gtf; allow overlap
bash run.saturation.20250823.2.sh

# 四、鉴定ORF
nohup bash call-orfs/merge.bam.20250923.sh &> ../log/merge.bam.20250923.log &
nohup bash call-orfs/ribo.seq.tools.build.anno.20250924.sh &> ../log/ribo.seq.tools.build.anno.20250924.log &
- 激活biotools conda环境之后PRICE在终端运行完毕
- 此外Warning: the stop codon is discontinuous,only first region is used, PB.17993.17这条warining可以忽略
nohup bash call-orfs/call.orfs.test.sh $(realpath ../processed/orfs) \
  $(realpath ../processed/merged_bam/Aligned.toTranscriptome.out.bam) \
  $(realpath ../processed/merged_bam/Aligned.sortedByCoord.out.bam)  &> ../log/call.orfs.20250923.log &

# 五、鉴定ORF之前transcriptome bam需要重新比对
nohup bash run.all.file.remap.newStarIndex.20250925.sh &> run.all.file.remap.newStarIndex.20250925.log &