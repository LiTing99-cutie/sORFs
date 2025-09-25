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
## 1.2 RNA-seq
## 1.3 检查比对结果
# 2.鉴定ORF
## 2.2 鉴定ORFs
## 2.3 检查鉴定结果
## 2.4.1 合并所有的bam [chosen]
# 2.2.Run.combine.call.orfs.v2.20250328.sh
bash 2.3.Uni.Organize_res_v2.20250331.sh human_brain_ribo_merge_call_orfs_20250338 0.05 1 custom
bash S3.1.Uni.Merge_Filter_Annotate.v1.20250331.sh

# 20250424 找出所有不管是经典还是非经典的小肽，方便后续做比较
bash 2.3.Uni.Organize_res_v3.20240424.sh human_brain_ribo_merge_call_orfs_20250338 0.05 1 custom 1 organized_include_cano
bash 3.1.Uni.Merge_Filter_Annotate.v2.20250424.sh $PWD/human_brain_ribo_merge_call_orfs_20250338/organized_include_cano/ \
    $PWD/human_brain_ribo_merge_call_orfs_20250338/merge_include_cano/
Rscript 3.2.rm_dup.filter.v2.20240425.R

# 20250520 更新参数
find $PWD/human_brain_output_20250227/ -name output | \
parallel -j 10 --joblog log/Run.call.orfs.change_ribocode_cutoff.log 'log_path={//}/log;bash 2.1.Uni.call.orfs.change_ribocode_cutoff.v1.20250521.sh {} Ribo_ORFs_change_ribocode_cutoff &> $log_path/Run.call.orfs.change_ribocode_cutoff.log'