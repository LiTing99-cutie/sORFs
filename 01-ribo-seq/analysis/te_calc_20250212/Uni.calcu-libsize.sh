#!/usr/bin/sh

################################################
#File Name: calcu-libsize.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Wed 31 Jul 2024 10:18:36 AM CST
################################################

# 该脚本的目的是从bam中得到reads的数量，如果是双端数据，只得到read 1的数量即可
# 20250221 增加sam的不同参数
set -eo pipefail

bam_lst=$1
output_path=$2
mkdir -p $output_path
if [ "$1" == "-h" ]; then
    echo "Usage: $0 <bam_lst> <output_path>"
    exit 0
fi
[ -f $output_path/libsize.txt ] || rm -rf $output_path/libsize.txt
for bam in $(cat $bam_lst);do
suffix="${bam##*.}"
if [ "$suffix" == "bam" ]; then
    lib=$(sambamba view -t 40 -c $bam)
else 
    lib=$(sambamba view -S -t 40 -c $bam)
fi
# 【修改点】匹配SRR或者ERR开头后面接数字，如果不符合，得不到对应结果
sample=$(basename $bam | grep -o -E 'ERR[0-9]+|SRR[0-9]+')
echo -e "$sample\t$lib" >> $output_path/libsize.txt
done