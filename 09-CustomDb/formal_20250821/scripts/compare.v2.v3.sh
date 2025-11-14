	#!/bin/bash

file1="/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103/compare_with_v2.orf_id.txt"
file2="/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/candidateORF.filtered.rmdup.renamed.pep.txt"

echo "=== 文件统计 ==="
echo ""

# 统计file1
if [ -f "$file1" ]; then
    count1=$(wc -l < "$file1")
    echo "v3 (filtered_v3_20251103): $count1 个ORF_id"
else
    echo "v3文件不存在: $file1"
    exit 1
fi

# 统计file2
if [ -f "$file2" ]; then
    count2=$(wc -l < "$file2")
    echo "v2 (filtered_v2_20251021): $count2 个ORF_id"
else
    echo "v2文件不存在: $file2"
    exit 1
fi

echo ""
echo "=== 重叠统计 ==="

# 统计重叠
overlap=$(comm -12 <(sort "$file1") <(sort "$file2") | wc -l)
echo "重叠 (两个文件都有): $overlap 个ORF_id"

# 统计v3独有
v3_only=$(comm -23 <(sort "$file1") <(sort "$file2") | wc -l)
echo "v3独有 (只在v3中): $v3_only 个ORF_id"

# 统计v2独有
v2_only=$(comm -13 <(sort "$file1") <(sort "$file2") | wc -l)
echo "v2独有 (只在v2中): $v2_only 个ORF_id"

echo ""
echo "=== 百分比 ==="
if [ $count1 -gt 0 ]; then
    overlap_pct1=$(echo "scale=2; $overlap * 100 / $count1" | bc)
    echo "v3中重叠占比: ${overlap_pct1}%"
fi

if [ $count2 -gt 0 ]; then
    overlap_pct2=$(echo "scale=2; $overlap * 100 / $count2" | bc)
    echo "v2中重叠占比: ${overlap_pct2}%"
fi

echo ""
echo "验证: $count1 (v3) = $overlap (重叠) + $v3_only (v3独有)"
echo "验证: $count2 (v2) = $overlap (重叠) + $v2_only (v2独有)"