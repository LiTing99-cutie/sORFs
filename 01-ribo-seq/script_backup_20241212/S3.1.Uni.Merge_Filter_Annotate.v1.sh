#!/usr/bin/sh
# 合并所有样本以及所有方法的结果

# 之前的结果输出文件夹名
orfs_res_dir=$1 # Ribo_ORFs_loose_para_20241206
# 合并后结果输出路径
output_path=$2 # ./mouse_brain_output_20241011/Ribo_ORFs_loose_para_20241206_merge

merge_res_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1a.Uni.merge_res.raw.v1.sh

mkdir -p $output_path
# 0. 加上样本名方法名，合并所有样本和方法的结果
find ./mouse_brain_output_20241011/ -name $orfs_res_dir | parallel -j 10 bash $merge_res_script {}
find ./mouse_brain_output_20241011/ -path "*/$orfs_res_dir/merge/nonCano.sorf.meta.merge.raw.3_ways.txt"|xargs cat > \
	$output_path/nonCano.sorf.meta.merge.raw.3_ways.all_samples.txt

# 1. 过滤：多个转录本编码同一个ORF，选择证据最强的转录本；去掉和已知注释蛋白序列相同的sORF；输出至少在一个样本和一个方法中鉴定得到的sORF
script_1=./Run.Select_repre_trans.R
Rscript $script_1
# 输出 nonCano.sorf.txt nonCano.sorf.filtered.txt

# 2. 增加元信息
script_2=./Run.add_meta.R
Rscript $script_2
# 输出 nonCano.sorf.filtered.add_meta.txt
script_3=./Run.annotate_ORF_type.R
Rscript $script_3
# 输出 nonCano.sorf.filtered.add_meta.orf_type_annotated.txt