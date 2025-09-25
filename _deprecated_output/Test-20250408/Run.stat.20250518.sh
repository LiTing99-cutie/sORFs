
#!/usr/bin/sh

################################################
#File Name: 1.3.Run.Mapping_Statistics.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 27 Mar 2025 04:13:27 PM CST
################################################

set -eo pipefail
output_dir=$PWD/01-output/stat
mkdir -p $output_dir

##### 0. 得到样本的不同输出list #####
find $PWD/ -name fq.stat.txt > $output_dir/fq.stat.lst
find $PWD/ -name *_Log.final.out > $output_dir/Log.final.out.lst
find $PWD -name "*_Aligned.sortedByCoord.out.bam" > $output_dir/Aligned.sortedByCoord.out.bam.lst

##### 1. 质控步骤过滤掉的reads百分比 #####
###### 1.0 原始READS数量 ######
ls /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/organize_all_test_data_20250515/*fq.gz |\
 egrep 0321 > $output_dir/raw_fastq.lst
cat $output_dir/raw_fastq.lst |xargs seqkit stat > $output_dir/raw_reads.txt
###### 1.1 被其他小RNA污染的比例 ######
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20241221-revision/Uni.Run.Mapping_Statistics.20250518.R
Rscript $script $output_dir/fq.stat.lst $output_dir

##### 2. 唯一比对数量以及唯一比对率 #####
source /home/user/data3/lit/project/sORFs/sORFs.utils.sh
extract_info $output_dir/Log.final.out.lst "Uniquely mapped reads %" \
    "$output_dir/Uniquely_mapped_reads_rate.txt"
extract_info $output_dir/Log.final.out.lst "Uniquely mapped reads number" \
    "$output_dir/Uniquely_mapped_reads_number.txt"
merge_files_by_sample "$output_dir/Uniquely_mapped_reads_rate.txt" "$output_dir/Uniquely_mapped_reads_number.txt" \
	"$output_dir/Uniquely_mapped_reads_rate_number.txt"

###### 3. 读段的长度分布 ######
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/1.3c.dataset.ribo.Rlenth.sh
nohup bash $script \
	$output_dir/Aligned.sortedByCoord.out.bam.lst \
	$output_dir/reads_length_distribution.txt &
###### 4. 读段在基因组上的分布 ######
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20241221-revision/Uni.read_genome_distri.20250518.sh
head -n1 $output_dir/Aligned.sortedByCoord.out.bam.lst > $output_dir/Aligned.sortedByCoord.out.bam.mouse.lst
tail -n1 $output_dir/Aligned.sortedByCoord.out.bam.lst > $output_dir/Aligned.sortedByCoord.out.bam.human.lst
bash $script $output_dir/Aligned.sortedByCoord.out.bam.mouse.lst $output_dir /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.bed
bash $script $output_dir/Aligned.sortedByCoord.out.bam.human.lst $output_dir
###### 5. 总体的三碱基周期性 ######
define_annotation_gencode_v41_human
conda activate ribocode
for bam in $(find $PWD -name "*Aligned.sortedByCoord.out.bam");do
[ -f $bam.bai ] || samtools index $bam
ribotish quality -b $bam -g $gtf -p 30
done

## 小鼠
define_annotation_gencode_v41_mouse
sample=E16_0321
bam=$PWD/output/${sample}.raw/output/alignment/${sample}.raw_Aligned.sortedByCoord.out.bam
ribotish quality -b $bam -g $gtf -p 30
# 编码基因所有codon的不同frame的比例
script=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/ribosome_profiling/parse_ribotish_qual.py
mkdir -p 01-output/stat/ribotish/
for sample in E16_0321.raw p21_0321.raw;do
echo "$sample"
python $script --sample_name $sample \
     --txt_path $PWD/output/$sample/output/alignment/${sample}_Aligned.sortedByCoord.out_qual.txt \
	 --offset_path $PWD/output/$sample/output/alignment/${sample}_Aligned.sortedByCoord.out.bam.para.py \
	 --RPF_start_distr_file $output_dir/ribotish/$sample.RPF_start_distr_file.txt \
	 --RPF_stop_distr_file $output_dir/ribotish/$sample.RPF_stop_distr_file.txt \
	 --frame_distr_file $output_dir/ribotish/$sample.frame_distr_file.txt
done
find ./ -name "*frame_distr_file.txt" |xargs cat > $output_dir/ribotish/Frame_distr_file.txt

###### 6. 整理和合并 ######
###### 24-36长度中三碱基周期性大于等于50%且显著性小于0.01的reads数量 ######
# 对于同一个样本，计算不同层次的Ribo-seq Reads结果
conda activate ribocode
# 得到24-36长度中三碱基周期性大于等于50%且显著性小于0.01的reads长度对应的bam
for sample in E16_0321.raw p21_0321.raw;do
pushd $PWD/output/$sample/output/Ribo-ORFs/RiboCode
bam_1=../../alignment/${sample}_Aligned.toTranscriptome.out.bam
bam_2=../../alignment/${sample}_Aligned.sortedByCoord.out.bam
if [ "$sample" = "E16_0321.raw" ]; then
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/RiboCode_annot/mm39
else
RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
fi
metaplots -a $RiboCode_annot -r $bam_1 -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o config_0.5_0.01
grep "^#" config_0.5_0.01_pre_config.txt|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.tab.txt
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250509/extract_specific_length_reads.sh
bash $script $bam_2 offset.tab.txt ../../alignment/${sample}_Aligned.sortedByCoord.out.withPeri.bam
popd
done
find ./ -name "*withPeri.bam" | xargs -I {} sh -c 'cnt=$(samtools view -c "{}"); echo -e "{}\t$cnt"' > $output_dir/p_sites_number.txt
###### 鉴定出的ORF数量 ######
find ./ -name "*collapsed.txt" |xargs wc -l  