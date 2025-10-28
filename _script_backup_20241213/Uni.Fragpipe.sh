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

/rd1/user/lit/software/fragpipe/bin/fragpipe --headless \
--workflow $workflow \
--manifest $manifest \
--workdir $workdir \
--ram 180 \
--threads $thread_n \
--config-tools-folder /rd1/user/lit/project/sORFs/fragpipe_tools \
--config-diann /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
--config-python ~/rd1/anaconda3/envs/py3/bin/python
