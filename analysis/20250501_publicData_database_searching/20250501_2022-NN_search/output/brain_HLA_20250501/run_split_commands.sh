#!/usr/bin/sh

################################################
#File Name: run_split_commands.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2025年05月06日 星期二 14时46分41秒
################################################

set -eo pipefail

# 原始命令模板（注意用双引号包裹变量）
SCRIPT_TEMPLATE='bash /rd1/user/lit/project/sORFs/Uni.Fragpipe.sh config/fragpipe.workflow part_file ./split_20_all_files_20250506 40 &> ./log/split_20_all_files_20250506/human.hla.cali.0.split20.part_suffix.log'

# 遍历所有拆分后的 part 文件
for part_file in config/fragpipe-files.fp-manifest.part-*; do
    # 提取 part 后缀（如 aa, ab）
    part_suffix=$(basename "$part_file" | cut -d'-' -f4)
    
    # 替换模板中的 {} 为实际 part 文件名和日志后缀
    cmd=$(echo "$SCRIPT_TEMPLATE" | sed "s|part_file|$part_file|g"|sed "s|part_suffix|$part_suffix|g")
    
    # 打印并执行命令（先打印验证，确认无误后去掉 echo 执行）
    echo "Running: $cmd"
    eval "$cmd"
done