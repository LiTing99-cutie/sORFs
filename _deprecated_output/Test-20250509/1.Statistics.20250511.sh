
#!/usr/bin/sh

################################################
#File Name: 1.3.Run.Mapping_Statistics.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 27 Mar 2025 04:13:27 PM CST
################################################

set -eo pipefail
output_dir=01-output/stat
mkdir -p $output_dir
##### 1. 质控步骤过滤掉的reads百分比 #####
ls $(find $PWD/ -name fq.stat.txt) \
    /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/output/fq.stat.txt > $output_dir/fq.stat.lst
Rscript 1a.Run.Mapping_Statistics.20250511.R

##### 2. 唯一比对数量以及唯一比对率 #####
ls $(find $PWD/ -name *_Log.final.out) \
	$(find /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/ -name *_Log.final.out) \
	> $output_dir/Log.final.out.lst
# 定义提取信息的函数
extract_info() {
    local search_lst=$1
    local search_term=$2
    local output_file=$3

    cat $search_lst| xargs grep "$search_term" |\
    awk -F'[:\t]' '{split($1, a, "/"); split(a[length(a)], b, "_Log.final.out"); print b[1] "\t" $NF}' > "$output_file"
}
merge_files_by_sample(){
	# 定义输入和输出文件路径
	rate_file=$1
	number_file=$2
	output_file=$3

	# 按第一列（样本名）合并，并添加表头
	echo -e "Sample\tUniquely_mapped_reads_rate\tUniquely_mapped_reads_number" > "$output_file"
	join -t $'\t' -1 1 -2 1 $rate_file $number_file >> "$output_file"
}
extract_info $output_dir/Log.final.out.lst "Uniquely mapped reads %" \
    "$output_dir/Uniquely_mapped_reads_rate.txt"
extract_info $output_dir/Log.final.out.lst "Uniquely mapped reads number" \
    "$output_dir/Uniquely_mapped_reads_number.txt"
merge_files_by_sample "$output_dir/Uniquely_mapped_reads_rate.txt" "$output_dir/Uniquely_mapped_reads_number.txt" \
	"$output_dir/Uniquely_mapped_reads_rate_number.txt"

###### 3. 读段的长度分布 ######
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/1.3c.dataset.ribo.Rlenth.sh

ls $(find $PWD -name "*_Aligned.sortedByCoord.out.bam") \
	$(find /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/ -name *_Aligned.sortedByCoord.out.bam) \
	> $output_dir/Aligned.sortedByCoord.out.bam.lst
nohup bash $script \
	$output_dir/Aligned.sortedByCoord.out.bam.lst \
	$output_dir/reads_length_distribution.txt &
###### 4. 读段在基因组上的分布 ######
bash 1b.Run.dataset.ribo.distri.20250511.sh
###### 5. 总体的三碱基周期性 ######
# 编码基因所有codon的不同frame的比例
script=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/ribosome_profiling/parse_ribotish_qual.py
mkdir -p 01-output/stat/ribotish/
for sample in p21_40_1_0425 p21_40_0425 p21_40_0422;do
echo "$sample"
python $script --txt_path $PWD/01-output/call-orfs/$sample.raw/output/alignment/$sample.raw_Aligned.sortedByCoord.out_qual.txt \
	 --offset_path $PWD/01-output/call-orfs/$sample.raw/output/alignment/$sample.raw_Aligned.sortedByCoord.out.bam.para.py \
	 --RPF_start_distr_file $PWD/01-output/stat/ribotish/$sample.RPF_start_distr_file.txt \
	 --RPF_stop_distr_file $PWD/01-output/stat/ribotish/$sample.RPF_stop_distr_file.txt \
	 --frame_distr_file $PWD/01-output/stat/ribotish/$sample.frame_distr_file.txt
done

sample=p21_0321
out_dir=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408
python $script --txt_path $out_dir/output/$sample.raw/output/alignment/$sample.raw_Aligned.sortedByCoord.out_qual.txt \
	 --offset_path $out_dir/output/$sample.raw/output/alignment/$sample.raw_Aligned.sortedByCoord.out.bam.para.py \
	 --RPF_start_distr_file $PWD/01-output/stat/ribotish/$sample.RPF_start_distr_file.txt \
	 --RPF_stop_distr_file $PWD/01-output/stat/ribotish/$sample.frame_distr_file.txt \
	 --frame_distr_file $PWD/01-output/stat/ribotish/$sample.frame_distr_file.txt

find ./ -name "p21*frame_distr_file.txt" |xargs cat > $PWD/01-output/stat/ribotish/total.frame_distr_file.txt
# 查看ribotish导出的pdf
ls $(find $PWD -name "*_Aligned.sortedByCoord.out_qual.pdf") \
	$(find /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/output/p21_0321.raw/ -name *_Aligned.sortedByCoord.out_qual.pdf)
###### 6. 整理和合并 ######

###### 原始READS数量 ######
ls /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250509/raw_data/*raw.fastq.gz /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/p21_0321.raw.fastq.gz > ./01-output/stat/raw_fastq.lst
cat ./01-output/stat/raw_fastq.lst |xargs seqkit stat > ./01-output/stat/raw_reads.txt

###### 24-36长度中三碱基周期性大于等于50%且显著性小于0.01的reads数量 ######
# 对于同一个样本，计算不同层次的Ribo-seq Reads结果
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
conda activate ribocode
# 得到24-36长度中三碱基周期性大于等于50%且显著性小于0.01的reads长度对应的bam
for sample in p21_40_1_0425 p21_40_0422 p21_40_0425;do
pushd 01-output/call-orfs/$sample.raw/output/Ribo-ORFs/RiboCode
metaplots -a $RiboCode_annot -r ../../alignment/$sample.raw_Aligned.toTranscriptome.out.bam -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o config_0.5_0.01
bam=../../alignment/$sample.raw_Aligned.sortedByCoord.out.bam
grep "^#" config_0.5_0.01_pre_config.txt|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.tab.txt
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250509/extract_specific_length_reads.sh
bash $script $bam offset.tab.txt ../../alignment/$sample.raw_Aligned.sortedByCoord.out.withPeri.bam
popd
done
find ./ -name "*withPeri.bam" | xargs -I {} sh -c 'cnt=$(samtools view -c "{}"); echo -e "{}\t$cnt"' > 01-output/stat/p_sites_number.txt

# 得到所有的比对（去掉比对不上的）
# 得到所有的指定长度的比对

###### 鉴定出的ORF数量 ######
find ./ -name "*collapsed.txt"