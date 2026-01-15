#!/usr/bin/sh

################################################
#File Name: Run.2024-11-01.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月01日 星期五 22时57分10秒
################################################

set -eo pipefail

source activate fragpipe

workflow=$1
manifest=$2
workdir=$3
thread_n=$4
bin_dir=/home/user/data2/lit/software/fragpipe/bin
tool_dir=/home/user/data2/lit/software/fragpipe_tools
$bin_dir/fragpipe --headless \
    --workflow $workflow \
    --manifest $manifest \
    --workdir $workdir \
    --ram 200 \
    --threads $thread_n \
    --config-tools-folder $tool_dir \
    --config-diann $tool_dir/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
    --config-python /home/user/data3/lit/anaconda3/envs/fragpipe/bin/python
