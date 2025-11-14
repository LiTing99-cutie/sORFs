#!/bin/bash

# ========================================
# ParaAT Ka/Ks 计算脚本 (修正版)
# ========================================

set -e
set -u
set -o pipefail

# ========================================
# 使用说明
# ========================================
usage() {
    cat << EOF
用法: $0 [选项]

必需参数:
  -p, --pep FILE        蛋白质序列文件 (.pep.fa)
  -c, --cds FILE        CDS序列文件 (.nucl.fa)
  -o, --output DIR      输出目录

可选参数:
  -m, --mode MODE       同源基因对生成模式 [human|all]
                        human: 第一个序列(人类) vs 其他所有 (默认)
                        all: 所有两两配对
  -a, --aligner TOOL    比对软件 [mafft|muscle|clustalw2] (默认: mafft)
  -t, --threads NUM     并行线程数 (默认: 10)
  -k, --kaks METHOD     Ka/Ks计算方法 [NG|LWL|MLWL|MA|MS|etc.] (默认: NG)
  -h, --help            显示此帮助信息

示例:
  $0 -p gene.pep.fa -c gene.nucl.fa -o output_dir
  $0 -p gene.pep.fa -c gene.nucl.fa -o output_dir -m all -k MA

EOF
    exit 1
}

# ========================================
# 默认参数
# ========================================
MODE="human"
ALIGNER="mafft"
THREADS=10
KAKS_METHOD="NG"
PEP=""
CDS=""
OUTPUT=""

# ========================================
# 解析命令行参数
# ========================================
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--pep)
            PEP="$2"
            shift 2
            ;;
        -c|--cds)
            CDS="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        -a|--aligner)
            ALIGNER="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -k|--kaks)
            KAKS_METHOD="$2"
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

# ========================================
# 检查必需参数
# ========================================
if [ -z "$PEP" ] || [ -z "$CDS" ] || [ -z "$OUTPUT" ]; then
    echo "错误: 缺少必需参数"
    usage
fi

# ========================================
# 检查文件存在性
# ========================================
if [ ! -f "$PEP" ]; then
    echo "错误: 蛋白质序列文件不存在: $PEP"
    exit 1
fi

if [ ! -f "$CDS" ]; then
    echo "错误: CDS序列文件不存在: $CDS"
    exit 1
fi

if [ "$MODE" != "human" ] && [ "$MODE" != "all" ]; then
    echo "错误: MODE参数必须是 'human' 或 'all'"
    exit 1
fi

# ========================================
# 创建输出目录
# ========================================
OUTPUT=$(realpath $OUTPUT)
mkdir -p $OUTPUT/paraAT_res

# ========================================
# 日志函数
# ========================================
LOG_FILE="${OUTPUT}/kaks_calculation.log"
[ -f $LOG_FILE ] && rm -rf $LOG_FILE
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a $LOG_FILE
}

# ========================================
# 开始处理
# ========================================
log "======================================"
log "Ka/Ks 计算流程开始"
log "======================================"
log "蛋白质序列: $PEP"
log "CDS序列: $CDS"
log "输出目录: $OUTPUT"
log "配对模式: $MODE"
log "比对软件: $ALIGNER"
log "线程数: $THREADS"
log "Ka/Ks方法: $KAKS_METHOD"
log ""

# ========================================
# 1. 检查序列数量
# ========================================
log "=== 步骤1: 检查序列 ==="
n_pep=$(grep -c "^>" $PEP)
n_cds=$(grep -c "^>" $CDS)

log "蛋白质序列数: $n_pep"
log "CDS序列数: $n_cds"

if [ $n_pep -ne $n_cds ]; then
    log "警告: PEP和CDS序列数量不一致！"
fi

if [ $n_pep -lt 2 ]; then
    log "错误: 需要至少2条序列才能进行配对"
    exit 1
fi

# ========================================
# 2. 生成同源基因对文件
# ========================================
log ""
log "=== 步骤2: 生成同源基因对 (模式: $MODE) ==="
HOMOLOG="${OUTPUT}/homolog.txt"

if [ "$MODE" == "human" ]; then
    # 第一遍找出hg38，第二遍生成配对
    seqkit seq -n $PEP | awk '
    /^hg38/ || /hg38/ {
        hg38 = $0
    }
    {
        seqs[NR] = $0
    }
    END {
        if (hg38 == "") {
            print "错误: 未找到hg38序列" > "/dev/stderr"
            exit 1
        }
        for (i = 1; i <= NR; i++) {
            if (seqs[i] != hg38) {
                print hg38 "\t" seqs[i]
            }
        }
    }' > $HOMOLOG
    n_pairs=$(wc -l < $HOMOLOG)
    log "生成 $n_pairs 对配对 (参考序列 vs 其他)"
    
elif [ "$MODE" == "all" ]; then
    seqkit seq -n $PEP | awk '{ids[NR]=$0} END {
        for(i=1; i<=NR-1; i++) {
            for(j=i+1; j<=NR; j++) {
                print ids[i]"\t"ids[j]
            }
        }
    }' > $HOMOLOG
    n_pairs=$(wc -l < $HOMOLOG)
    log "生成 $n_pairs 对配对 (所有两两组合)"
fi

log "配对列表前5行:"
head -5 $HOMOLOG | tee -a $LOG_FILE

# ========================================
# 3. 创建进程控制文件
# ========================================
log ""
log "=== 步骤3: 配置并行处理 ==="
PROC_FILE="$OUTPUT/proc"
echo $THREADS > $PROC_FILE
log "使用 $THREADS 个线程"

# ========================================
# 4. 运行ParaAT
# ========================================
log ""
log "=== 步骤4: 运行ParaAT ==="

[ -d $OUTPUT/paraAT_res ] && rm -rf $OUTPUT/paraAT_res
mkdir -p $OUTPUT/paraAT_res

log "ParaAT参数:"
log "  比对软件: $ALIGNER"
log "  输出格式: axt"
log "  移除gap: 是"
log "  移除错配: 是"
log "  计算KaKs: 是"
[ -d $OUTPUT/paraAT_res ] && rm -rf $OUTPUT/paraAT_res
ParaAT.pl \
    -h $HOMOLOG \
    -n $CDS \
    -a $PEP \
    -processor $PROC_FILE \
    -m $ALIGNER \
    -f axt \
    -g \
    -t \
    -o $OUTPUT/paraAT_res 2>&1 | tee -a $LOG_FILE

if [ $? -ne 0 ]; then
    log "错误: ParaAT运行失败"
    exit 1
fi

log "ParaAT完成"

# ========================================
# 5. 计算Ka/Ks
# ========================================
log ""
log "=== 步骤5: 计算 Ka/Ks (方法: $KAKS_METHOD) ==="

axt_count=$(find $OUTPUT/paraAT_res -name "*.axt" | wc -l)
log "找到 $axt_count 个axt文件"

if [ $axt_count -eq 0 ]; then
    log "错误: 未找到axt文件"
    exit 1
fi

for axt in $OUTPUT/paraAT_res/*.axt; do
    base=$(basename $axt .axt)
    log "处理: $base"
    
    KaKs_Calculator -i $axt -o $OUTPUT/paraAT_res/${base}.kaks -m $KAKS_METHOD 2>&1 | tee -a $LOG_FILE
    
    if [ $? -ne 0 ]; then
        log "警告: $base 计算失败"
    fi
done

log "Ka/Ks计算完成"

# ========================================
# 6. 汇总结果 - 详细版（使用实际的列）
# ========================================
log ""
log "=== 步骤6: 汇总结果 ==="

KAKS_RES="${OUTPUT}/kaks.res"
# 详细结果 - 保留所有原始列并添加选择类型
log "生成详细结果: $KAKS_RES"

# 写入表头
echo -e "Species_pair\tKa\tKs\tKa/Ks\tS-Substitutions\tN-Substitutions" > $KAKS_RES

for kaks in $OUTPUT/paraAT_res/*aln.kaks; do
    tail -n +2 $kaks | awk '{
        print $1"\t"$3"\t"$4"\t"$5"\t"$12"\t"$13
    }'
done >> $KAKS_RES