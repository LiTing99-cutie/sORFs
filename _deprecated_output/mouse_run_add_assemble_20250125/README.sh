#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 25 Jan 2025 01:03:15 PM CST
################################################

set -eo pipefail

# 复制后根据修改gtf和genepred
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/S2.2.Uni.loose.parameter.sh Uni.run.tree_tools.sh

# 复制后修改genepred
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.0.Uni.Organize_res_v2.sh S3.0.Uni.Organize_res_v2.sh

cp /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1.Uni.Merge_Filter_Annotate.v1.sh S3.1.Uni.Merge_Filter_Annotate.v1.sh
# 运行
mkdir -p log
nohup bash Run.20250125.add_assemble.sh &> log/Run.20250125.add_assemble.log &