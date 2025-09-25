#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 27 Feb 2025 03:47:18 PM CST
################################################

set -eo pipefail

# 这个路径的目的是重复在鼠上面的分析，分为几个步骤
# 1、比对，包括RNA-seq和Ribo-seq质控后fastq的比对，以及基于RNA-seq的转录本的组装
# 2、组装完成后，基于三个软件的ORF的鉴定
# 3、合并和过滤

mkdir log
# 1、比对（Ribo-seq还需要进一步去除污染）
## 1.1 Ribo-seq
### 基本质控后的fastq路径
find /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized -name "*fq.gz" |grep qc/ > trimmed.ribo.fq.lst
for fq in $(cat trimmed.ribo.fq.lst);do
    echo "***Processing $fq"
	bash 1.1.Uni.ribo.rmContam.mapping.human.20250227.sh $fq $PWD/human_brain_output_20250227
done &> log/ribo.rmContam.mapping.human.20250227.1.log

# 内存不够，终端了，重新跑
tail -n +27 trimmed.ribo.fq.lst > trimmed.ribo.fq.remaining.lst
for fq in $(cat trimmed.ribo.fq.remaining.lst);do
    echo "***Processing $fq"
	bash 1.1.Uni.ribo.rmContam.mapping.human.20250227.sh $fq $PWD/human_brain_output_20250227
done &> log/ribo.rmContam.mapping.human.20250227.2.log

## 1.2 RNA-seq
# 由于当时基本质控的时候并不是同一批跑的，所以目录结构并没有统一，因此这里我们把输出统一到一个文件夹内
cat <(find /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized -name "*fq.gz" |grep qc_rna_seq/) \
	<(find /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/2022-NN/qc -name "*fq.gz") \
	<(find /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized -path "*/RNA/*fastq.gz" |egrep "2022-MC|2014-JN") > trimmed.rna.fq.lst
nohup bash 1.2.Uni.rna.mapping.assemble.human.20250227.sh \
	$PWD/trimmed.rna.fq.lst \
	$PWD/human_brain_rna_seq_alignment_assemble &> log/ribo.rmContam.mapping.human.20250227.log &
nohup bash 1.2.1.Run.rna.mapping.assemble.human.20250228.padding.sh &

## 1.3 检查比对结果
bash 1.3.Run.Mapping_Statistics.sh

# 2.鉴定ORF
## 2.1 将组装的转录本和基因加到参考的基因注释中 【这一步没有做，根据在鼠中的分析，引入这一步虽然能有一些额外的结果，但是最终呈现方式并没有想好，且可能会对结果有影响】
nohup bash 2.0.Run.Build_annotation.sh &> log/Run.Build_annotation.log &
## 2.2 鉴定ORFs
find $PWD/human_brain_output_20250227/ -name output | \
parallel -j 10 --joblog log/Run.call.orfs.log 'log_path={//}/log;bash 2.1.Uni.call.orfs.sh {} Ribo_ORFs &> $log_path/Run.call.orfs.log'

## 2.3 检查鉴定结果
### 查看有多少样本是无法鉴定出p-site位置的；26个样本
find ./ -name Run.call.orfs.log |xargs grep -H "Error, can not determine" |grep -oP "SRR\d+" > human_brain_output_20250227/2022_NN_ribocode_cannot_determine.srr.txt
### 查看有多少样本退出状态为0；28个样本
less log/Run.call.orfs.log | tail -n +2|awk '$7!=0'|grep -oP "SRR\d+"|uniq > human_brain_output_20250227/failed.srr.txt
### overlap为24个
cat human_brain_output_20250227/2022_NN_ribocode_cannot_determine.srr.txt human_brain_output_20250227/failed.srr.txt|sort|uniq -c|awk '$1==2'|wc -l
### 退出状态为0的还有4个是因为内存不足等问题
grep -v -f human_brain_output_20250227/2022_NN_ribocode_cannot_determine.srr.txt human_brain_output_20250227/failed.srr.txt
# SRR15513149
# SRR15906442
# SRR15906444
# SRR15906424
## 这个是因为memory不足
find ./ -name Run.call.orfs.log |grep -E "SRR15513149|"|xargs grep -H memory
# 这个可能是因为并行运行的问题
find ./ -name Run.call.orfs.log |grep -E "SRR15906444|SRR15906424|SRR15906442"|xargs grep java.lang.RuntimeException

grep -v -f human_brain_output_20250227/failed.srr.txt human_brain_output_20250227/2022_NN_ribocode_cannot_determine.srr.txt
# SRR15906448
# SRR15906464

# 失败了，但是并不是因为ribocode无法发现周期性的拿去重跑
grep -v -f human_brain_output_20250227/2022_NN_ribocode_cannot_determine.srr.txt human_brain_output_20250227/failed.srr.txt > human_brain_output_20250227/consider_rerun.txt

ls -ld human_brain_output_20250227/*/* | cut -f3 -d "/"|head -n 84 > human_brain_output_20250227/all_run.txt
# 从所有的样本中去掉单独无法call ORF的样本
grep -v -f human_brain_output_20250227/2022_NN_ribocode_cannot_determine.srr.txt human_brain_output_20250227/all_run.txt > human_brain_output_20250227/high_quality_run.txt
find ./ -name metaplots_pre_config.txt |grep -f human_brain_output_20250227/high_quality_run.txt|xargs cat > human_brain_output_20250227/high_quality_run.meta.config.txt

## 2.4 重新run之前失败的四个样本；以及合并所有的bam文件，进行合并运行
nohup bash Run.combine.call.orfs.20250312.sh &> log/Run.combine.call.orfs.20250312.log &
nohup bash Run.combine.call.orfs.20250312.sh &> log/Run.combine.call.orfs.20250312.1.log &

### 重新run之前失败的四个样本
find $PWD/human_brain_output_20250227/ -name output | grep -f human_brain_output_20250227/consider_rerun.txt | \
	parallel -j 3 --joblog log/Run.call.orfs.rerun.20250327.log 'log_path={//}/log;bash 2.1.Uni.call.orfs.sh {} Ribo_ORFs &> $log_path/Run.call.orfs.rerun.20250327.log'
# SRR15906442运行PRICE步骤特别慢
find $PWD/human_brain_output_20250227/ -name output | grep SRR15906442 | \
	parallel --tmpdir /home/user/data3/lit/tmp -j 1 --joblog log/Run.call.orfs.rerun.20250328.log 'log_path={//}/log;bash 2.1.Uni.call.orfs.sh {} Ribo_ORFs &> $log_path/Run.call.orfs.rerun.20250328.log'
# 这里面只有SRR15513149重跑成功
find ./ -name Run.call.orfs.rerun.20250327.log |grep -E "SRR15906444|SRR15906424|SRR15513149|SRR15906442"

echo "SRR15513149" | cat - human_brain_output_20250227/high_quality_run.txt > human_brain_output_20250227/high_quality_run.add.20250328.txt
egrep -v "SRR15906444|SRR15906424|SRR15906442" human_brain_output_20250227/high_quality_run.add.20250328.txt > human_brain_output_20250227/high_quality_run.add.rm.20250328.txt

# 20250328
# 其实有些样本的结果并不包含在最后的结果log/Run.call.orfs.log中，因此，我们只需要找到ribocode鉴定有周期性的特定读长的reads就可以了。

# Run.combine.call.orfs.20250312.sh重命名为Run.combine.call.orfs.v1.20250312.sh

## 2.4.1 合并所有的bam [chosen]
# 2.2.Run.combine.call.orfs.v2.20250328.sh

# 由于RibORF没有跑完，先整理其他两个软件的结果
bash 2.3.Uni.Organize_res_v2.20250331.sh human_brain_ribo_merge_call_orfs_20250338 0.05 1 custom
/home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1.Uni.Merge_Filter_Annotate.v1.20250331.sh