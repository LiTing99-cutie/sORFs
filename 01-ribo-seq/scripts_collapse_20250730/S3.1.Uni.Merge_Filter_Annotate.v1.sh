#!/usr/bin/sh
# 合并所有方法的结果


# PRICE
awk -v OFS='\t' '{print $1,$2,"PRICE"}' PRICE/nonCano.sorf.tab > nonCano.sorf.meta.raw.PRICE.txt

# RiboCode
awk -v OFS='\t' '{print $1,$2,"RiboCode"}' RiboCode/nonCano.sorf.tab > nonCano.sorf.meta.raw.RiboCode.txt

# RibORF
awk -v OFS='\t' '{print $1,$2,"RibORF"}' RibORF/nonCano.sorf.tab > nonCano.sorf.meta.raw.RibORF.txt

# merge
cat nonCano.sorf.meta.raw.PRICE.txt nonCano.sorf.meta.raw.RiboCode.txt nonCano.sorf.meta.raw.RibORF.txt | \
awk -v OFS='\t' '{print $1,$2,$3,"'$sample'"}' > nonCano.sorf.meta.merge.raw.3_ways.txt

# 之前的结果输出文件夹名
orfs_res_dir=$1 # Ribo_ORFs_loose_para_20241206
# 合并后结果输出路径
output_path=$2 # ./mouse_brain_output_20241011/Ribo_ORFs_loose_para_20241206
proj_path=/home/user/data3/lit/project/sORFs/01-ribo-seq
merge_res_script=$proj_path/S3.1a.Uni.merge_res.raw.v1.sh
script_1=$proj_path/S3.1b.Run.Select_repre_trans.v1.1.R
script_2=$proj_path/S3.1c.Run.add_meta.v1.1.R
script_3=$proj_path/S3.1d.Run.annotate_ORF_type.v1.1.R

mkdir -p $output_path/merge $output_path/filter
# 0. 加上样本名方法名，合并所有样本和方法的结果
find $proj_path/mouse_brain_output_20241011/ -path "*/output/$orfs_res_dir" | parallel -j 10 bash $merge_res_script {}
find $proj_path/mouse_brain_output_20241011/ -path "*/$orfs_res_dir/merge/nonCano.sorf.meta.merge.raw.3_ways.txt"|xargs cat > \
	$output_path/merge/nonCano.sorf.meta.merge.raw.3_ways.all_samples.txt

# 1. 过滤：多个转录本编码同一个ORF，选择证据最强的转录本；去掉和已知注释蛋白序列相同的sORF；输出至少在一个样本和一个方法中鉴定得到的sORF
Rscript $script_1 $output_path/merge/nonCano.sorf.meta.merge.raw.3_ways.all_samples.txt $output_path/filter

# 2. 增加元信息
Rscript $script_2 $output_path/filter/nonCano.sorf.filtered.txt $output_path/filter/nonCano.sorf.meta.merge.raw.3_ways.all_samples.topTrans.txt $output_path/filter/nonCano.sorf.filtered.add_meta.txt

Rscript $script_3 $output_path/filter/nonCano.sorf.filtered.add_meta.txt $output_path/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt 