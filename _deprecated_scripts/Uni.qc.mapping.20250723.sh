# 参数
raw_fastq_lst=$1
output_path=$2
# 脚本路径定义
script_1=Uni.qc.single.v1.20250723.sh
script_2=Uni.qc.batch.v1.20250723.sh
script_3=Uni.rmContam.mapping.v1.20250723.sh
# 创建输出文件夹
mkdir -p $output_path/fastqc $output_path/alignment

# 1. qc
bash $script_2 $output_path/fastqc \
<(cat $raw_fastq_lst) \
yes yes yes $script_1 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
# 2. remove contaminant and mapping
for trimmed_fq in $(find $output_path -name "*trimmed.fq.gz");do
bash $script_3 $trimmed_fq $output_path/alignment
done