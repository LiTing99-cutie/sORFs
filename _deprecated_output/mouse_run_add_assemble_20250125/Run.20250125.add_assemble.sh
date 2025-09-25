#!/usr/bin/sh

################################################
#File Name: Run.20250125.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Sat 25 Jan 2025 04:02:58 PM CST
################################################

set -eo pipefail

mkdir -p log
proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/
# export到环境变量方便parallel运行
export run_name=Ribo_ORFs_add_assemble_20250125
export date=20250125
export script_1=S2.Uni.run.tree_tools.sh
export script_2=S3.0.Uni.Organize_res_v2.sh
export script_3=S3.1.Uni.Merge_Filter_Annotate.v1.sh
export script_4=$proj_path/S3.2.Visual.1.1.R
export pre_output_path=$proj_path/mouse_brain_output_20241011
uniprot_fasta=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform_2024_10_25.fasta
# build annotation
echo -e "build annotation at $(date '+%Y-%m-%d %H:%M:%S')"
# already done
# bash S1.Build_annotation.sh
# call ORFs
echo -e "call ORFs at $(date '+%Y-%m-%d %H:%M:%S')"
# find $pre_output_path -name output | \
# parallel -j 10 --joblog log/Run.tree_tools.$date.log 'log_path={//}/log;bash $script_1 {} $run_name &> $log_path/Run.tree_tools.$date.log'
# organize results
echo -e "organize results at $(date '+%Y-%m-%d %H:%M:%S')"
# find $pre_output_path -path "*/output/$run_name" | \
# parallel -j 10 --joblog log/Organize_res.$date.log 'bash $script_2 {} 0.05 0 custom'
# merge,filter and add meta columns
echo -e "merge,filter and add meta columns at $(date '+%Y-%m-%d %H:%M:%S')"
bash $script_3 $run_name $pre_output_path/$run_name
# visualization
echo -e "visualization at $(date '+%Y-%m-%d %H:%M:%S')"
output_path=$pre_output_path/$run_name
mkdir -p $output_path/visual/
Rscript $script_4 $output_path/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt \
	$output_path/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.filtered.txt  $output_path/visual/

# 生成fa文件去搜库
echo -e "Generate database at $(date '+%Y-%m-%d %H:%M:%S')"
mkdir -p $proj_path/database/custom/$run_name && pushd $proj_path/database/custom/$run_name
output_path=$pre_output_path/$run_name
less $output_path/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt | tail -n +2 | awk -v OFS='\t' '{print $2,$4}' |seqkit tab2fx > nonCano.sorf.$date.fa
cat $uniprot_fasta nonCano.sorf.$date.fa > uniprot.nonCano.sorf.$date.fa
popd