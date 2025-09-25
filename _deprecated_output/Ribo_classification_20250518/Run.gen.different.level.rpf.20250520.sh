mkdir -p output/bam/bam_1
mkdir -p output/bam/bam_2
mkdir -p output/fq
# 合并所有的特异性比对的bam
bam_lst=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/Aligned.sortedByCoord.out.bam.lst
bam_trans_lst=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/Aligned.toTranscriptome.out.bam.lst
extract_specific_length_reads_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/in_house_phase_I_data_20250520/Uni.extract_specific_length_reads.v1.20250520.sh
samtools merge -@ 40 -f output/bam/bam_1/merged.bam -b $bam_lst
samtools merge -@ 40 -f output/bam/bam_1/merged.trans.bam -b $bam_trans_lst
samtools view -@ 40 -h output/bam/bam_1/merged.bam > output/bam/bam_1/merged.sam
# 选择24-36的片段
bash $extract_specific_length_reads_script output/bam/bam_1/merged.bam output/bam/bam_2/length.24-36.bam
# 所有长度在24-36且三碱基周期性大于等于50%的核糖体保护片段【P-site】
get_p_sites_v1(){
    source activate ribocode
    # 输入的bam文件，需要绝对路径
    bam=$1
    # 输入的sam文件
    sam_input=$2
    # 输出路径
    output_path=$3
    mkdir -p $output_path
    # 输出文件名
    sam_output=$4
    peri_cutoff=$5
    p_value_cutoff=$6
    RiboCode_annot=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/annotation/RiboCode_annot
    RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
    # 进入输出路径
    pushd $output_path
    metaplots -a $RiboCode_annot -r $bam -f0_percent $peri_cutoff -pv1 $p_value_cutoff -pv2 $p_value_cutoff
    grep "^#" metaplots_pre_config.txt|awk -v OFS='\t' '{print $2,$4+3}'|tail -n +3|head -n -1 > offset.correction.parameters.txt
    perl $RibORF_path/offsetCorrect.pl -r $sam_input -p offset.correction.parameters.txt -o $sam_output
    # 转换为bam，更加节省空间
    bam_output=$(echo $sam_output|sed 's/.sam/.bam/')
    samtools view -@ 40 -bh $sam_output > $bam_output && rm -rf $sam_output
    popd
}
export -f get_p_sites_v1
## 生成参数列表
: > para.lst
for bam_trans in $(cat $bam_trans_lst);do
sample=$(basename -s "_Aligned.toTranscriptome.out.bam" $bam_trans)
sam=$(echo $bam_trans|sed 's/_Aligned.toTranscriptome.out.bam/.sam/;s/data3/data/')
output_dir=$PWD/output/bam/bam_3/$sample
output_file=$sample.p_sites_0.5_1.sam
echo "$bam_trans $sam $output_dir $output_file"
done >> para.lst
## 使用 parallel 并行运行
parallel --joblog log/get_p_sites_0.5_1.log -j 10 --colsep ' ' get_p_sites_v1 {1} {2} {3} {4} 0.5 1 :::: para.lst
samtools merge -@ 40 -f output/bam/bam_3/merged.bam $(find output/bam/bam_3/ -name "*p_sites_0.5_1.bam")

# 所有长度在24-36且三碱基周期性大于等于50%且第一个阅读框的P site数目显著大于第二个和第三个阅读框的核糖体保护片段【P-site】
: > para.1.lst
for bam_trans in $(cat $bam_trans_lst);do
sample=$(basename -s "_Aligned.toTranscriptome.out.bam" $bam_trans)
sam=$(echo $bam_trans|sed 's/_Aligned.toTranscriptome.out.bam/.sam/;s/data3/data/')
output_dir=$PWD/output/bam/bam_4/$sample
output_file=$sample.p_sites_0.5_0.01.sam
echo "$bam_trans $sam $output_dir $output_file" >> para.1.lst
done 
parallel -j 10 --colsep ' ' get_p_sites_v1 {1} {2} {3} {4} 0.5 0.01 :::: para.1.lst
samtools merge -@ 40 -f output/bam/bam_4/merge.bam $(find output/bam/bam_4/ -name "*p_sites_0.5_0.01.bam")


