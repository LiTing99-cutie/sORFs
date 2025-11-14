#!/bin/bash
#
# calculate_rpkm_tpm.sh (优化版)
# 从RibORF结果计算RPKM、TPM和密码子覆盖度
#
# 用法: 
#   bash calculate_rpkm_tpm.sh -i input.txt -o output.txt

set -euo pipefail

# ============ 参数解析 ============
usage() {
    cat << EOF
用法: $0 -i <input_file> -o <output_file> [-h]

必需参数:
  -i, --input     输入文件 (RibORF结果文件)
  -o, --output    输出文件

可选参数:
  -h, --help      显示帮助信息

输入文件格式:
  - TSV格式，可以有#开头的注释行
  - 列1: ORF_id
  - 列6: ORF长度(bp)
  - 列7: RPF reads数
  - 列8: P-sites数

输出文件包含:
  ORF_id, RPF_reads, Psites_number, RPF_RPKM, Psites_RPKM, 
  RPF_TPM, Psites_TPM, RPF_codon_coverage, Psites_codon_coverage

示例:
  $0 -i riborf_result.txt -o rpkm_tpm_output.txt

EOF
    exit 1
}

# 初始化变量
INPUT_FILE=""
OUTPUT_FILE=""

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "错误: 未知参数 $1"
            usage
            ;;
    esac
done

# 检查必需参数
if [[ -z "$INPUT_FILE" ]] || [[ -z "$OUTPUT_FILE" ]]; then
    echo "错误: 缺少必需参数"
    usage
fi

# 检查输入文件是否存在
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "错误: 输入文件不存在: $INPUT_FILE"
    exit 1
fi

# 检查输出目录是否存在，不存在则创建
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# ============ 主程序 ============
echo "===== 计算RPKM、TPM和密码子覆盖度 (优化版) ====="
echo "输入文件: $INPUT_FILE"
echo "输出文件: $OUTPUT_FILE"
echo ""

# 计算总reads数（用于RPKM计算）
echo "步骤1: 计算总reads数..."
read rpf_total ps_total < <(
    grep -v '^#' "$INPUT_FILE" | awk 'NR==1{next} {rpf+=$7; ps+=$8} END{print rpf, ps}'
)
: "${rpf_total:=0}"
: "${ps_total:=0}"

echo "  RPF总reads数: $rpf_total"
echo "  P-sites总数: $ps_total"

# 检查是否有数据
if [[ "$rpf_total" -eq 0 ]] && [[ "$ps_total" -eq 0 ]]; then
    echo "警告: 未检测到有效数据"
fi

# 第一遍扫描：计算RPK总和
echo ""
echo "步骤2: 计算RPK总和..."
read rpf_rpk_sum ps_rpk_sum < <(
    grep -v '^#' "$INPUT_FILE" | awk '
    NR==1{next}
    {
        len_bp=$6+0; rpf=$7+0; ps=$8+0;
        if(len_bp>0){
            rpf_rpk_sum += (rpf*1000)/len_bp;
            ps_rpk_sum  += (ps*1000)/len_bp;
        }
    }
    END{print rpf_rpk_sum, ps_rpk_sum}'
)
: "${rpf_rpk_sum:=0}"
: "${ps_rpk_sum:=0}"

echo "  RPF_RPK总和: $rpf_rpk_sum"
echo "  Psites_RPK总和: $ps_rpk_sum"

# 第二遍扫描：计算所有指标并输出
echo ""
echo "步骤3: 计算RPKM、TPM和密码子覆盖度..."
grep -v '^#' "$INPUT_FILE" | awk -v rpfT="$rpf_total" -v psT="$ps_total" \
    -v rpf_rpk_sum="$rpf_rpk_sum" -v ps_rpk_sum="$ps_rpk_sum" '
BEGIN{
    OFS="\t";
    # 输出表头
    print "ORF_id","RPF_reads","Psites_number","RPF_RPKM","Psites_RPKM","RPF_TPM","Psites_TPM","RPF_codon_coverage","Psites_codon_coverage"
}
NR==1{next}
{
    id=$1; len_bp=$6+0; rpf=$7+0; ps=$8+0;
    
    # 跳过空行或无效行
    if(id=="" || len_bp<=0) next;
    
    # RPKM计算: (reads × 10^9) / (length_bp × total_reads)
    rpf_rpkm = (rpfT>0) ? (rpf*1e9)/(len_bp*rpfT) : 0;
    ps_rpkm  = (psT>0)  ? (ps *1e9)/(len_bp*psT)  : 0;
    
    # TPM计算: (RPK / sum(RPK)) × 10^6
    rpf_rpk = (rpf*1000)/len_bp;
    ps_rpk  = (ps*1000)/len_bp;
    rpf_tpm = (rpf_rpk_sum>0) ? (rpf_rpk/rpf_rpk_sum)*1e6 : 0;
    ps_tpm  = (ps_rpk_sum>0)  ? (ps_rpk/ps_rpk_sum)*1e6   : 0;
    
    # 每密码子覆盖度（len_bp为bp，1个密码子=3bp）
    rpf_cov = (rpf*3)/len_bp;
    ps_cov  = (ps*3)/len_bp;
    
    print id, rpf, ps, rpf_rpkm, ps_rpkm, rpf_tpm, ps_tpm, rpf_cov, ps_cov;
    
    count++;
}
END{
    # 输出统计信息到stderr
    print "处理了 "count" 个ORF" > "/dev/stderr"
}' > "$OUTPUT_FILE"

# 检查输出
if [[ -f "$OUTPUT_FILE" ]]; then
    line_count=$(wc -l < "$OUTPUT_FILE")
    echo ""
    echo "===== 完成 ====="
    echo "输出文件: $OUTPUT_FILE"
    echo "输出行数: $line_count (包含表头)"
    echo ""
    echo "列说明:"
    echo "  1. ORF_id              - ORF标识符"
    echo "  2. RPF_reads           - RPF reads数"
    echo "  3. Psites_number       - P-sites数"
    echo "  4. RPF_RPKM            - RPF的RPKM值"
    echo "  5. Psites_RPKM         - P-sites的RPKM值"
    echo "  6. RPF_TPM             - RPF的TPM值"
    echo "  7. Psites_TPM          - P-sites的TPM值"
    echo "  8. RPF_codon_coverage  - RPF每密码子覆盖度"
    echo "  9. Psites_codon_coverage - P-sites每密码子覆盖度"
    echo ""
    echo "前5行预览:"
    head -n 5 "$OUTPUT_FILE" | column -t -s $'\t'
else
    echo "错误: 输出文件创建失败"
    exit 1
fi