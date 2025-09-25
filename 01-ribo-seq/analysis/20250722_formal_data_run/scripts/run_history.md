# qc，去除污染，比对
## 生成原始数据list文件
raw_data_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/in_house_human_organized_20250625
ls $raw_data_path/*.fq.gz > ../processed/raw_fastq_batch_1.20250723.lst
nohup bash Uni.qc.mapping.20250723.sh ../processed/raw_fastq_batch_1.20250723.lst ../processed/batch_1 &> ../log/Run.qc.mapping.20250723.log &
- 备注：不小心把p21_0523_{1..4}也给跑了，所以等到跑完之后再删除吧
- 第二个步骤失败，因为没有输入绝对路径

nohup bash Uni.qc.mapping.20250723.1.sh /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/raw_fastq_batch_1.20250723.lst /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1 &> ../log/Run.qc.mapping.20250723.1.log &
- 使用绝对路径进行输入
- bowtie比对比较慢，15个样本预计要2.2天
- 这部分结果比较慢且中间因为服务器内存不够有些样本结果不完整，已经删除

nohup bash star.rmContam.20250724.sh &> ../log/STAR.rm.rRNA.20250724.log &
- 测试使用STAR替换
- 该样本同样使用bowtie2的50个线程需要运行1小时28分钟，现在测试同样使用STAR需要多久，以及结果是否相同
- 根据../processed/batch_1/alignment/p21_0626_2/output/fq.stat.txt，bowtie2过滤之后还剩下294,696,873条reads
- 294,696,873/303,960,902
- 287,676,696/303,960,902
- 287,675,694 (减少seedSearchStartLmax参数)
- 287,966,415 (不允许软剪切以及splice-aware比对)

trimmed_fastq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1/fastqc/p21_0626_2/output/trimmed_fastq/p21_0626_2_trimmed.fq.gz
nohup bash Uni.rmContam.mapping.using.star.v2.20250724.sh $trimmed_fastq ../tmp/test_for_p21_0626_2 &> ../log/test_for_p21_0626_2.20250724.log &
- 进一步测试使用STAR替换，优化参数

nohup bash Uni.qc.mapping.using.star.v2.20250724.sh /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/raw_fastq_batch_1.20250723.lst /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1 &> ../log/Run.qc.mapping.using.star.v2.20250724.log &
- 使用STAR替换
- 由于STAR比对的参数有所更新，且index也有所更新，这部分结果也删除

nohup bash star.rmContam.20250725.sh &> ../log/STAR.rm.rRNA.20250725.log &
- 分别使用新的rRNA index文件（U替换成T）以及RNA central+新的rRNA index文件（U替换成T）的rRNA构建的序列

nohup bash Uni.qc.mapping.using.star.v3.20250725.sh /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/raw_fastq_batch_1.20250723.lst /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1 &> ../log/Run.qc.mapping.using.star.v3.20250725.log &
- 使用STAR替换，但是更换了rRNA的index以及outFilterMultimapNmax参数

nohup bash run.extract.25-34.length.reads.20250725.sh &
- 提取特异性比对的reads中符合特定长度的reads

nohup bash Uni.saturation.calcu.gene.counts.v1.20250728.sh ../processed/batch_1/filtered_bam ../processed/batch_1/geneCounts/ &> ../log/saturation.calcu.gene.counts.20250728.log &
- gene counts以及transcript counts

bash Uni.transcript_coverage_analysis_v1_20250725.sh ../processed/batch_1/geneCounts/transcript.counts.txt  ../results/batch_1/coverage_analysis/transcript
bash Uni.transcript_coverage_analysis_v1_20250725.sh ../processed/batch_1/geneCounts/gene.counts.txt  ../results/batch_1/coverage_analysis/gene
- 进行饱和度分析

bash Uni.transcript_coverage_analysis_v1_20250725.sh ../processed/batch_1/geneCounts/gene.counts.txt  ../results/batch_1/coverage_analysis/gene_permu_1

# 评估ribo-seq文库结果
bash ribo.qual.assess.20250730.sh
Rscript Uni.organize.v1.20250730.R /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/results/batch_1/qual_assess \
    /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/figures/batch_1/qual_assess


