#!/usr/bin/env bash
set -euo pipefail

# ==================== 配置部分 ====================
# 脚本路径配置
script_dir="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis"
Mapping_Statistics_script="$script_dir/Test-20241221-revision/Uni.Run.Mapping_Statistics.20250518.R"
calcu_ribo_len_script="$script_dir/Run_for_human_20250227/1.3c.dataset.ribo.Rlenth.sh"
read_genome_distri_script="$script_dir/Test-20241221-revision/Uni.read_genome_distri.20250518.sh"
parse_ribotish_qual_script="/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/ribosome_profiling/parse_ribotish_qual.py"
RiboCode_annot="$script_dir/Run_for_human_20250227/annotation/RiboCode_annot"
extract_specific_length_reads_script="$script_dir/Test-20250509/extract_specific_length_reads.sh"

# 输出目录配置
output_dir="$PWD/01-output/stat"
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
find "$PWD" -name "fq.stat.txt" > "$output_dir/fq.stat.lst"
find "$PWD" -name "*_Log.final.out" > "$output_dir/Log.final.out.lst"
find "$PWD" -name "*_Aligned.sortedByCoord.out.bam" > "$output_dir/Aligned.sortedByCoord.out.bam.lst"

# 1. 质控统计
echo "处理质控数据..."
# 1.0 原始reads数量
check_command seqkit
cat "$PWD/raw_fastq.lst" | xargs seqkit stat > "$output_dir/raw_reads.txt"

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
bash "$read_genome_distri_script" "$output_dir/Aligned.sortedByCoord.out.bam.lst" "$output_dir"

# 5. 三碱基周期性分析
echo "分析三碱基周期性..."
source activate ribocode
check_command ribotish
check_command samtools

# 为所有BAM文件创建索引
for bam in $(find "$PWD" -name "*Aligned.sortedByCoord.out.bam"); do
    [ -f "$bam.bai" ] || samtools index "$bam"
    ribotish quality -b "$bam" -g "$gtf" -p 30
done

# 解析ribotish结果
for sample in $(ls "$PWD/01-output/call-orfs/"); do
    echo "处理样本: $sample"
    input_txt="$PWD/01-output/call-orfs/$sample/output/alignment/${sample}_Aligned.sortedByCoord.out_qual.txt"
    input_offset="$PWD/01-output/call-orfs/$sample/output/alignment/${sample}_Aligned.sortedByCoord.out.bam.para.py"
    
    check_file "$input_txt"
    check_file "$input_offset"
    
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

# 6. 高级分析
echo "进行高级分析..."

# 6.1 提取特定长度读段
for sample in $(ls "$PWD/01-output/call-orfs/"); do
    echo "处理样本: $sample"
    pushd "$PWD/01-output/call-orfs/$sample/output/Ribo-ORFs/RiboCode" > /dev/null
    
    bam_1="../../alignment/${sample}_Aligned.toTranscriptome.out.bam"
    bam_2="../../alignment/${sample}_Aligned.sortedByCoord.out.bam"
    
    check_file "$bam_1"
    check_file "$bam_2"
    
    metaplots -a "$RiboCode_annot" -r "$bam_1" -f0_percent 0.5 -pv1 0.01 -pv2 0.01 -o config_0.5_0.01
    grep "^#" config_0.5_0.01_pre_config.txt | awk -v OFS='\t' '{print $2,$4+3}' | tail -n +3 | head -n -1 > offset.tab.txt
    
    bash "$extract_specific_length_reads_script" "$bam_2" offset.tab.txt "../../alignment/${sample}_Aligned.sortedByCoord.out.withPeri.bam"
    
    popd > /dev/null
done

# 6.2 统计p位点数量
find "$PWD" -name "*withPeri.bam" | xargs -I {} sh -c 'cnt=$(samtools view -c "{}"); echo -e "{}\t$cnt"' > "$output_dir/p_sites_number.txt"

# 6.3 统计ORF数量
find "$PWD" -name "*collapsed.txt" | while read -r file; do
    lines=$(wc -l < "$file")
    name=$(basename "$file" "_collapsed.txt")
    printf "%s\t%s\n" "$name" "$lines"
done > "$output_dir/orf_number.txt"

echo "====== 处理完成 ======"
echo "结果保存在: $output_dir"