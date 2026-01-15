#!/bin/bash
source activate ribocode

# ================= 路径与参数配置 =================
# 输入 BAM (绝对路径，这点很好)
INPUT_BAM="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed/merged_bam/Aligned.toTranscriptome.out.bam"

# 项目路径
PROJ_PATH="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf"
RIB_ANNOT="${PROJ_PATH}/processed/annotation/RiboCode_annot"

# 输出主目录 (为了防止 cd 后相对路径失效，这里获取脚本运行时的绝对路径)
# 如果 ../processed 存在，这样写最稳妥
OUT_DIR="$(readlink -f ../processed/saturation_analysis)" 
# 或者如果你的环境不支持 readlink -f，就写死绝对路径，比如 /home/user/.../processed/saturation_analysis
mkdir -p $OUT_DIR

SEED=2025
f0_percent_cutoff=0.6
PVALUE_CUTOFF=0.01
# ===============================================

echo "开始执行饱和度分析 (10% -> 100%)..."

for i in {1..10}; do
    PERCENT=$((i * 10))
    echo "========================================"
    echo "正在处理测序深度: ${PERCENT}%"
    
    # 1. 创建子目录
    SAMPLE_DIR="${OUT_DIR}/depth_${PERCENT}"
    mkdir -p $SAMPLE_DIR
    
    # === 关键修正 A：进入子目录 ===
    # 这样 metaplots 生成的所有文件都会乖乖待在 depth_xx 里面
    cd $SAMPLE_DIR || exit 1
    
    # 注意：进入子目录后，文件名只需写当前目录下的名字
    SUBSAMPLE_BAM="subset_${PERCENT}.bam"
    
    # 2. 准备 BAM 文件
    if [ $i -eq 10 ]; then
        echo "深度 100%: 链接原始 BAM 文件..."
        ln -sf $INPUT_BAM $SUBSAMPLE_BAM
    else
        echo "深度 ${PERCENT}%: 执行 samtools 抽样..."
        # 40个线程对于view可能有点多(瓶颈在IO)，但如果服务器强劲也没问题
        samtools view -@ 40 -b -s ${SEED}.${i} $INPUT_BAM > $SUBSAMPLE_BAM
    fi

    # 3. 运行 Metaplots
    echo "Running metaplots..."
    
    # 不需要 -o，因为已经在子目录里了
    metaplots -a $RIB_ANNOT \
              -r $SUBSAMPLE_BAM \
              -f0_percent $f0_percent_cutoff \
              -pv1 $PVALUE_CUTOFF \
              -pv2 $PVALUE_CUTOFF 

    # 4. 运行 RiboCode
    
    # === 关键修正 B：定义变量 ===
    # metaplots 默认生成的固定文件名
    CONFIG_FILE="metaplots_pre_config.txt"
    
    if [ -f "$CONFIG_FILE" ]; then
        echo "Running RiboCode..."
        RiboCode --alt_start_codons CTG,ACG,GTG,TTG \
                 --pval-cutoff 0.05 \
                 -a $RIB_ANNOT \
                 -c $CONFIG_FILE \
                 -l no -g \
                 -o "RiboCode_result" \
                 --output-gtf \
                 --output-bed \
                 --min-AA-length 6
                 
        echo "深度 ${PERCENT}% 分析完成。"
    else
        echo "错误: 未找到配置文件 $CONFIG_FILE，跳过此步骤。"
    fi
    
    # === 关键修正 C：退回上一级，准备下一次循环 ===
    cd - > /dev/null

done

echo "所有饱和度分析任务已完成。"