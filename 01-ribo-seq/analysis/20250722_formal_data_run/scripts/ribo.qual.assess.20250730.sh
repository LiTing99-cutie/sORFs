#!/usr/bin/env bash
set -euo pipefail

# ==================== 配置部分 ====================
# 脚本路径配置
script_dir="stat_scripts"
Mapping_Statistics_script="$script_dir/Uni.Run.Mapping_Statistics.20250518.R"
calcu_ribo_len_script="$script_dir/Calcu.dataset.ribo.Rlenth.sh"
read_genome_distri_script="$script_dir/Uni.read_genome_distri.20250518.sh"
parse_ribotish_qual_script="$script_dir/parse_ribotish_qual.py"
extract_specific_length_reads_script="Uni.extract_specific_length_reads.v1.20250520.sh"
parse_offdict_script="$script_dir/parse_offdict.py"
# 输出目录配置
# res_dir=$1
# 参数一：结果路径
res_dir=$(realpath ../processed/batch_1/)
bam_dir=$res_dir/filtered_bam
# 参数二：输出结果路径
output_dir=$(realpath ../results/batch_1/qual_assess)
# 参数三：原始fastq路径
raw_fastq_lst=../processed/raw_fastq_batch_1.20250723.lst
mkdir -p "$output_dir"
mkdir -p "$output_dir/ribotish"

# 工具配置
source "/home/user/data2/lit/bin/lit_utils.sh"
define_annotation_gencode_v41_human

# ==================== 功能函数 ====================
# 检查命令是否存在
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "错误: 未找到命令 $1"
        exit 1
    fi
}

# 检查文件是否存在
check_file() {
    if [ ! -f "$1" ]; then
        echo "错误: 文件 $1 不存在"
        exit 1
    fi
}

# ==================== 主流程 ====================
echo "====== 开始处理 ======"

# 0. 生成样本列表
echo "生成样本列表..."
find "$res_dir" -name "fq.stat.txt" > "$output_dir/fq.stat.lst"
find "$res_dir" -name "*_Log.final.out" > "$output_dir/Log.final.out.lst"
find "$bam_dir" -name "*_Aligned.sortedByCoord.out.bam" > "$output_dir/Aligned.sortedByCoord.out.bam.lst"

# 1. 质控统计
echo "处理质控数据..."
# 1.0 原始reads数量
check_command seqkit
cat "$raw_fastq_lst" | xargs seqkit stat -j 30 > "$output_dir/raw_reads.txt"

# 1.1 被其他小RNA污染的比例
check_command Rscript
check_file "$Mapping_Statistics_script"
Rscript "$Mapping_Statistics_script" "$output_dir/fq.stat.lst" "$output_dir"

# 2. 唯一比对统计
echo "处理比对数据..."
extract_info "$output_dir/Log.final.out.lst" "Uniquely mapped reads %" \
    "$output_dir/Uniquely_mapped_reads_rate.txt"
extract_info "$output_dir/Log.final.out.lst" "Uniquely mapped reads number" \
    "$output_dir/Uniquely_mapped_reads_number.txt"
merge_files_by_sample "$output_dir/Uniquely_mapped_reads_rate.txt" \
    "$output_dir/Uniquely_mapped_reads_number.txt" \
    "$output_dir/Uniquely_mapped_reads_rate_number.txt"

# 3. 读段长度分布
echo "分析读段长度分布..."
check_file "$calcu_ribo_len_script"
bash "$calcu_ribo_len_script" \
    "$output_dir/Aligned.sortedByCoord.out.bam.lst" \
    "$output_dir/reads_length_distribution.txt"

# 4. 基因组分布
echo "分析读段基因组分布..."
check_file "$read_genome_distri_script"
bash "$read_genome_distri_script" "$output_dir/Aligned.sortedByCoord.out.bam.lst" "$output_dir/r_distri"

# 5. 三碱基周期性分析
echo "分析三碱基周期性..."
source activate ribocode
check_command ribotish
check_command samtools

# ribotish
for bam in $(cat "$output_dir/Aligned.sortedByCoord.out.bam.lst"); do
    sample=$(basename $bam | sed 's/_Aligned.sortedByCoord.out.bam//')
    echo $sample
    [ -f "$bam.bai" ] || samtools index -@ 30 "$bam"
    cutoff=0.6
    ribotish quality -b "$bam" -g "$gtf" -p 30 --th 0.6 -r "$output_dir/ribotish/${sample}_0.6.para.py" \
        -f "$output_dir/ribotish/${sample}_0.6_qual.pdf" \
        -o "$output_dir/ribotish/${sample}_0.6_qual.txt"
    python $parse_offdict_script $output_dir/ribotish/${sample}_0.6.para.py > $output_dir/ribotish/${sample}_0.6_offset_tab.txt
done

# 解析ribotish结果
for file in $(ls $output_dir/ribotish/*_0.6_qual.txt); do
    sample=$(basename $file | sed 's/_0.6_qual.txt//')
    echo "处理样本: $sample"
    input_txt="$output_dir/ribotish/${sample}_0.6_qual.txt"
    input_offset="$output_dir/ribotish/${sample}_0.6.para.py"
    
    python "$parse_ribotish_qual_script" \
        --sample_name "$sample" \
        --txt_path "$input_txt" \
        --offset_path "$input_offset" \
        --RPF_start_distr_file "$output_dir/ribotish/$sample.RPF_start_distr_file.txt" \
        --RPF_stop_distr_file "$output_dir/ribotish/$sample.RPF_stop_distr_file.txt" \
        --frame_distr_file "$output_dir/ribotish/$sample.frame_distr_file.txt"
done

# 合并frame分布结果
find "$output_dir/ribotish/" -name "*frame_distr_file.txt" | xargs cat > "$output_dir/ribotish/Frame_distr_file.txt"

# p sites的数量
>"$output_dir"/ribotish/offset.tab.all.txt
for file in "$output_dir"/ribotish/*_0.6_offset_tab.txt; do
    sample=$(basename "$file" _0.6_offset_tab.txt)
    awk -v sample="$sample" -v OFS='\t' '{print sample, $1, $2}' "$file"
done |cat >> "$output_dir"/ribotish/offset.tab.all.txt

echo "====== 处理完成 ======"
echo "结果保存在: $output_dir"