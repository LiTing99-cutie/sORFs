#!/bin/bash
# 脚本名称：extract_rna_stats_with_path.sh
# 功能：提取样本信息并包含原始BAM路径

input_file=$1
output_file="enhanced_results_with_path.tsv"

echo -e "BamPath\tSampleID\tFailedFraction\tPlusPlus\tPlusMinus\tDifference\tLibraryType" > "$output_file"

awk -v OFS="\t" '
/\/[^/]+_Aligned\.sortedByCoord\.out\.bam$/ {
    # 保留原始路径
    bam_path = $0;
    
    # 提取样本ID
    sample_id = $0;
    sub(/.*\//, "", sample_id);
    sub(/_Aligned\.sortedByCoord\.out\.bam$/, "", sample_id);
    
    # 读取后续行
    for (i=1; i<=6; i++) {
        getline;
        if (i == 4) { failed = $NF }
        if (i == 5) { plusplus = $NF }
        if (i == 6) { plusminus = $NF }
    }
    
    # 计算和判断
    diff = plusminus - plusplus;
    if (diff > 0.2) {
        lib_type = 2;
    } else if (diff < -0.2) {
        lib_type = 1;
    } else {
        lib_type = 0;
    }
    
    # 输出结果（增加bamPath作为第一列）
    print bam_path, sample_id, failed, plusplus, plusminus, diff, lib_type;
}
' "$input_file" >> "$output_file"

echo "分析完成！带路径的结果已保存到 $output_file"