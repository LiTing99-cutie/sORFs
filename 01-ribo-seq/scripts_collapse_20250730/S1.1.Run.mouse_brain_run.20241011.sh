#!/usr/bin/sh

################################################
#File Name: mouse_brain_run.20241004.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 04 Oct 2024 11:56:50 AM CST
################################################

set -eo pipefail
# 20241011
#### 修改引用的脚本 ####

script=Ribo.Run.Uni.20241011.sh

for fq in /home/user/data3/licq/peptidomics/PublicData/mouse_brain/2020_NAR_WangH_E15.5_P42/ribo/*fastq.gz;do
	bash $script $fq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
done

for fq in /home/user/data3/licq/peptidomics/PublicData/mouse_brain/2015_Science_ChoJ_hippocampus/ribo/*fastq.gz;do
	bash $script $fq TGGAATTCTCGGGTGCCAAGG
done

for fq in /home/user/data3/licq/peptidomics/PublicData/mouse_brain/2019_multispecies_mouse_brain/ribo/*fastq.gz;do
	bash $script $fq AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
done

echo "All done"