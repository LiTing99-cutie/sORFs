#!/usr/bin/sh

################################################
#File Name: Run.2024-11-01.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月01日 星期五 22时57分10秒
################################################

set -eo pipefail

source activate fragpipe

mkdir -p output/MS/Fragpipe_output/2024-11-01-2 output/MS/Fragpipe_output/2024-11-01-3

/rd1/user/lit/software/fragpipe/bin/fragpipe --headless \
--workflow /rd1/user/lit/project/sORFs/config/2024-11-01-2-null-uniprot-custom/fragpipe.workflow \
--manifest /rd1/user/lit/project/sORFs/config/2024-11-01-2-null-uniprot-custom/fragpipe-files.fp-manifest \
--workdir /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024-11-01-2 \
--ram 100 \
--threads 30 \
--config-tools-folder /rd1/user/lit/project/sORFs/fragpipe_tools \
--config-diann /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
--config-python ~/rd1/anaconda3/envs/py3/bin/python

/rd1/user/lit/software/fragpipe/bin/fragpipe --headless \
--workflow /rd1/user/lit/project/sORFs/config/2024-11-01-3-trypsin-lysc-uniprot-custom/fragpipe.workflow \
--manifest /rd1/user/lit/project/sORFs/config/2024-11-01-3-trypsin-lysc-uniprot-custom/fragpipe-files.fp-manifest \
--workdir /rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024-11-01-03 \
--ram 100 \
--threads 30 \
--config-tools-folder /rd1/user/lit/project/sORFs/fragpipe_tools \
--config-diann /rd1/user/lit/project/sORFs/fragpipe_tools/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 \
--config-python ~/rd1/anaconda3/envs/py3/bin/python
