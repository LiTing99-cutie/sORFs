##### 重新比对，允许多重比对 #####

mkdir output log
trimmed_clean_fq_dir=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227
find $trimmed_clean_fq_dir -name "*trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz" >  trimmed_clean_fq.lst
fq_lst=$PWD/trimmed_clean_fq.lst
# script_dir=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227
for fq in $(cat $fq_lst);do
    echo "***Processing $fq"
	bash Uni.get_all_ribo_mapping_reads.20250518.sh $fq $PWD/output/all_ribo_mapping_reads_20250518
done &> log/all_ribo_mapping_reads_20250518.log

# 内存不够，剩下的再运行
# 从第 34 行开始处理
sed -n '34,$p' "$fq_lst" | while read fq; do
    echo "***Processing $fq"
    # 运行脚本并记录日志
    bash "Uni.get_all_ribo_mapping_reads.20250518.sh" "$fq" $PWD/output/all_ribo_mapping_reads_20250518 
done &> log/all_ribo_mapping_reads_20250520.log

sed -n '33p' "$fq_lst" | while read fq; do
    echo "***Processing $fq"
    # 运行脚本并记录日志
    bash "Uni.get_all_ribo_mapping_reads.20250518.sh" "$fq" $PWD/output/all_ribo_mapping_reads_20250518 
done &> log/all_ribo_mapping_reads_20250521.log

##### 降低三碱基周期性的阈值，鉴定ORF #####
# 每个样本单独鉴定
# see /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227
# 合并鉴定
bam_dir=$PWD/output/bam/bam_1
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/in_house_phase_I_data_20250520/Uni.combine.call.orfs.v1.20250520.sh
nohup bash $script $PWD/output/orfs $bam_dir/merged.trans.bam $bam_dir/merged.bam $bam_dir/merged.sam &> log/combine.call.orfs.log &
# PRICE的内存不够，重新运行
bam_dir=$PWD/output/bam/bam_1
script_rerun=Tmp.uni.combined.call.orfs.v1.20250520.sh
nohup bash $script_rerun $PWD/output/orfs $bam_dir/merged.trans.bam $bam_dir/merged.bam $bam_dir/merged.sam &> log/combine.call.orfs.rerun.20250522.log &
# 联系占用内存多的清理之后再运行
nohup bash $script_rerun $PWD/output/orfs $bam_dir/merged.trans.bam $bam_dir/merged.bam $bam_dir/merged.sam &> log/combine.call.orfs.rerun.20250523.log &
##### 生成level 1的RPF #####
mkdir output/bam/bam_0/
find ./ -name "*_Aligned.sortedByCoord.out.bam" > allow_mul_map.bam.lst
path=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Ribo_classification_20250518/output/bam/bam_0
mkdir -p $path
# 可能是内存不太够了，所以得到的merged.bam和所有的bam的大小不同
nohup samtools merge -@ 40 -f $path/merged.bam -b allow_mul_map.bam.lst &> log/merge.mul.bam.log &
##### 生成level 2到level 5的RPF #####
bash Run.gen.different.level.rpf.20250520.sh
##### 统计 #####
samtools view -@ 40 -c output/bam/bam_0/merged.bam 
# 576072955 576M
samtools view -@ 40 -c output/bam/bam_1/merged.bam 
# 498060943 498M 
samtools view -@ 40 -c output/bam/bam_2/length.24-36.bam 
# 393812469 393M
samtools view -@ 40 -c output/bam/bam_3/merged.bam
# 361109415 361M
samtools view -@ 40 -c output/bam/bam_4/merged.bam

# 把所有的/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Ribo_classification_20250518/output/all_ribo_mapping_reads_20250518/移动到相应的data目录下面，然后删除源目录下对应的文件