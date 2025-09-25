# 关键输入：去除污染之后的reads list；比对到基因组和转录组上的bam
# 对于同一个样本，可以合并bam进行call orfs的操作
mkdir -p output/bam/bam_1
mkdir -p output/bam/bam_2
mkdir -p output/fq
mkdir log
# 首先整理p21_0321，p21_40_0425，p21_40_1_0425，p21_40_0422这四个文库的去除污染之后的reads
pattern="*trimmed.rRNA.tRNA.snoRNA.unaligned.fq.gz"
(find $path_1 -name $pattern|grep p21_0321 ; find $path_2 -name $pattern) > output/fq/fq.lst
zcat $(cat output/fq/fq.lst) |gzip > output/fq/merged.fastq.gz
# 比对到基因组上，允许多重比对
cleaned_fastq=output/fq/merged.fastq.gz
STARindex=/home/user/data3/lit/resource/genome/human/hg38/hg38_STARindex_v2.5.2b/
sample=human_brain
conda activate biotools
[ -d output/alignment ] || mkdir output/alignment
STAR \
--readFilesIn $cleaned_fastq \
--readFilesCommand zcat \
--seedSearchStartLmax 15 \
--runThreadN 40 \
--outFilterMismatchNmax 2 \
--genomeDir $STARindex \
--outFileNamePrefix output/alignment/${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--outFilterMultimapNmax 9999 \
--outFilterMatchNmin 16 \
--outSAMattributes All \
--alignEndsType EndToEnd &> log/STAR.log
# 合并所有的特异性比对的bam
path_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408
path_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250509
pattern="*_Aligned.sortedByCoord.out.bam"
(find $path_1 -name $pattern|grep p21_0321 ; find $path_2 -name $pattern) > output/bam/bam_1/bam.lst
samtools merge -@ 40 -f output/bam/bam_1/merged.bam -b output/bam/bam_1/bam.lst
pattern="*_Aligned.toTranscriptome.out.bam"
(find $path_1 -name $pattern|grep p21_0321 ; find $path_2 -name $pattern) > output/bam/bam_1/bam_trans.lst
samtools merge -@ 40 -f output/bam/bam_1/merged.trans.bam -b output/bam/bam_1/bam_trans.lst
samtools view -@ 40 -h output/bam/bam_1/merged.bam > output/bam/bam_1/merged.sam
# 选择24-36的片段
bash Uni.extract_specific_length_reads.v1.20250520.sh output/bam/bam_1/merged.bam output/bam/bam_2/length.24-36.bam
# 使用三个工具来鉴定ORF
bam_dir=$PWD/output/bam/bam_1
bash Uni.combine.call.orfs.v1.20250520.sh $PWD/output/orfs $bam_dir/merged.trans.bam $bam_dir//merged.bam $bam_dir/merged.sam 

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
    popd
}
get_p_sites_v1 $PWD/output/bam/bam_1/merged.trans.bam \
    $PWD/output/bam/bam_1/merged.sam \
    $PWD/output/bam/bam_3/ \
    p_sites_0.5_1.sam \
    0.5 \
    1
# 所有长度在24-36且三碱基周期性大于等于50%且第一个阅读框的P site数目显著大于第二个和第三个阅读框的核糖体保护片段【P-site】
get_p_sites_v1 $PWD/output/bam/bam_1/merged.trans.bam \
    $PWD/output/bam/bam_1/merged.sam \
    $PWD/output/bam/bam_4/ \
    p_sites_0.5_0.01.sam \
    0.5 \
    0.01

# 1 output/alignment/human_brain_Aligned.sortedByCoord.out.bam
# 2 output/bam/bam_1/merged.bam
# 3 output/bam/bam_2/length.24-36.bam
# 4 output/bam/bam_3/p_sites_0.5_1.sam
# 5 output/bam/bam_4/p_sites_0.5_0.01.sam

samtools view -@ 40 -c output/alignment/human_brain_Aligned.sortedByCoord.out.bam 
samtools view -@ 40 -c output/bam/bam_1/merged.bam 
samtools view -@ 40 -c output/bam/bam_2/length.24-36.bam 
wc -l output/bam/bam_3/p_sites_0.5_1.sam 
wc -l output/bam/bam_4/p_sites_0.5_0.01.sam

# 计算指定的ORF的不同level的PRF的数量
bash Run.i.sorf.o.gpe_gtf_ribo_expr.v1.20250521.sh