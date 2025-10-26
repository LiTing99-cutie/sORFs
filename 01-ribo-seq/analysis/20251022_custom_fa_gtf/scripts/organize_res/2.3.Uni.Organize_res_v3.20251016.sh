#!/usr/bin/env bash
# S3.0.Uni.Merge_res_v1.sh — unify outputs from PRICE / RiboCode / RibORF
# Usage:
#   bash S3.0.Uni.Merge_res_v1.sh <work_path> <pvalue_cutoff> <only_cano_Scodon:0|1> <riborf_cutoff:custom|youden|fixed> <only_cano_type:0|1> <sub_dir> [min_aa] [max_aa]
#
# Example:
#   bash S3.0.Uni.Merge_res_v1.sh mouse_brain_output_20241011/*/*/output/Ribo*ORFs* 0.05 0 custom 0 merged 6 NA

set -euo pipefail
IFS=$'\n\t'

############################
# ------- 参数解析 -------- #
############################
work_path=${1:? "ERROR: missing <work_path>"}
pvalue_cutoff=${2:? "ERROR: missing <pvalue_cutoff> (e.g., 0.05)"}
only_cano_Scodon=${3:? "ERROR: missing <only_cano_Scodon> (0|1)"}
riborf_cutoff=${4:? "ERROR: missing <riborf_cutoff> (custom|youden|fixed)"}
only_cano_type=${5:? "ERROR: missing <only_cano_type> (0|1)"}
sub_dir=${6:? "ERROR: missing <sub_dir>"}

# 可选长度参数（默认不限制上限；min_aa 默认 6）
min_aa=${7:-6}
max_aa=${8:-NA}

############################
# ------- 常量路径 -------- #
############################
genePred=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.15.gpe
fa=/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/custom_ref.fa

project_path=/home/user/data3/lit/project/sORFs/01-ribo-seq
organize_res_script_path=$project_path/analysis/20250813_demo_data_analysis/scripts/organize_res
translate_script=$project_path/scripts_collapse_20250730/S3.0c.Uni.translate_gtf.v2.20250325.sh

price_script=$organize_res_script_path/S3.0b.Uni.Generate_genepred_PRICE_v1.py
compare_script_path=$project_path/ref/Ribo-seq-Tool-Comparison-Scripts-v2.0/Scripts_for_RiboCode_Analysis
riborf_script=$organize_res_script_path/S3.0a.Uni.Format_RibORF_v3.202501016.R
ribocode_script=$organize_res_script_path/S3.0f.Uni.Filter_RiboCode.20251016.R

sample=Human_brain  # 20250331 固定样本名（如需动态可在此处修改）

############################
# ------- 小工具函数 ------- #
############################
log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "command not found: $1"; }

pick_latest_file() {
  # 在通配符匹配的多个文件中选择修改时间最新的一个
  local pattern="$1"
  shopt -s nullglob
  local files=( $pattern )
  shopt -u nullglob
  ((${#files[@]})) || die "no file matched: $pattern"
  # 按 mtime 倒序取第一个
  ls -1t "${files[@]}" | head -n1
}

############################
# ------- 依赖检查 -------- #
############################
for c in python Rscript seqkit genePredToGtf samtools; do need "$c"; done
[[ -d "$work_path" ]] || die "work_path not found: $work_path"
[[ -s "$genePred" ]] || die "genePred not found: $genePred"
[[ -s "$fa" ]] || die "FASTA not found: $fa"
[[ -s "$price_script" ]] || die "price_script not found: $price_script"
[[ -s "$translate_script" ]] || die "translate_script not found: $translate_script"
[[ -s "$riborf_script" ]] || die "riborf_script not found: $riborf_script"
[[ -s "$ribocode_script" ]] || die "ribocode_script not found: $ribocode_script"

# Conda 环境（可选）
if command -v conda >/dev/null 2>&1; then
  # 使用 base 里的 python（与原脚本一致）
  # shellcheck disable=SC1091
  source activate base || true
fi

############################
# ------- 目录准备 -------- #
############################
cd "$work_path"
mkdir -p "$sub_dir"/{RibORF,PRICE,RiboCode}
log "Working under: $work_path"
log "Sub-dir: $sub_dir"
log "Params: pvalue_cutoff=$pvalue_cutoff only_cano_Scodon=$only_cano_Scodon riborf_cutoff=$riborf_cutoff only_cano_type=$only_cano_type min_aa=$min_aa max_aa=$max_aa"

################################
# ------- PRICE 模块 ---------- #
################################
log "== PRICE =="
cd "$sub_dir/PRICE"

# PRICE 输出：挑选最新的 *.orfs.tsv
price_orfs_tsv=$(pick_latest_file "../../PRICE/*.orfs.tsv")
log "PRICE orfs TSV: $price_orfs_tsv"

# 1) PRICE -> sORF genePred
python "$price_script" "$price_orfs_tsv" "$genePred" "$pvalue_cutoff" "$only_cano_Scodon" "$only_cano_type" nonCano.formatted.gpe
[[ -s nonCano.formatted.gpe ]] || die "PRICE formatted gpe not generated."

# 2) genePred -> GTF -> FASTA (蛋白)
genePredToGtf file nonCano.formatted.gpe nonCano.formatted.gtf
bash "$translate_script" nonCano.formatted.gtf "$fa" $PWD
mv -f prot.fa nonCano.fa
[[ -s nonCano.fa ]] || die "nonCano.fa not generated."

# 3) 按 AA 长度过滤到 sORF
if [[ -z "${max_aa}" || "${max_aa}" == "NA" ]]; then
  max_for_seqkit=-1
else
  max_for_seqkit="${max_aa}"
fi
seqkit seq -g -m "$min_aa" -M "$max_for_seqkit" nonCano.fa > nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
[[ -s nonCano.sorf.fa && -s nonCano.sorf.tab ]] || die "PRICE sORF fasta/tab not generated."

################################
# ------- RiboCode 模块 ------- #
################################
log "== RiboCode =="
cd ../RiboCode

# RiboCode 的主结果（样本同名）
ribocode_txt="../../RiboCode/${sample}.txt"
[[ -s "$ribocode_txt" ]] || die "RiboCode result not found: $ribocode_txt"

# 是否启用长度过滤：按照你的要求，这里**显式启用**并传入 min/max（max=NA 时在 R 脚本里不会限制上限）
# 若你想改成“默认不启用”，把下面的 1 改为 0 即可。
enable_len_filter=1
log "Run RiboCode filter: enable_len_filter=$enable_len_filter (min_aa=$min_aa, max_aa=$max_aa)"
Rscript "$ribocode_script" "$ribocode_txt" "$pvalue_cutoff" "$only_cano_Scodon" "$only_cano_type" nonCano.sorf.txt "$enable_len_filter"  "$min_aa" "$max_aa" 
[[ -s nonCano.sorf.txt ]] || die "RiboCode filtered table not generated."

# 生成 fasta/tab
python "$compare_script_path/Generate_Fasta_RiboCode.v1.20251016.py" nonCano.sorf.txt nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
[[ -s nonCano.sorf.fa && -s nonCano.sorf.tab ]] || die "RiboCode sORF fasta/tab not generated."

################################
# -------- RibORF 模块 -------- #
################################
log "== RibORF =="
cd ../RibORF

riborf_param="../../RibORF/repre.valid.pred.pvalue.parameters.txt"
riborf_gpe="../../RibORF/repre.valid.ORF.genepred.txt"
riborf_stat="../../RibORF/stat.cutoff.txt"
[[ -s "$riborf_param" ]] || die "RibORF param not found: $riborf_param"
[[ -s "$riborf_gpe"   ]] || die "RibORF genePred not found: $riborf_gpe"
[[ -s "$riborf_stat"  ]] || die "RibORF cutoff stat not found: $riborf_stat"

# 长度过滤：你的 R 脚本设计为“默认不启用；仅当提供 min_aa/max_aa 时启用”
# 此处我们**始终传入**两个参数（min_aa、max_aa），其中 max_aa=NA 时等价于“仅下限”
log "Run RibORF formatter with length bounds (min_aa=$min_aa, max_aa=$max_aa)"
Rscript "$riborf_script" "$riborf_param" "$riborf_gpe" "$riborf_stat" "$riborf_cutoff" "$only_cano_Scodon" "$only_cano_type" nonCano.sorf.formatted.gpe "$min_aa" "$max_aa" 
[[ -s nonCano.sorf.formatted.gpe ]] || die "RibORF formatted gpe not generated."

# genePred -> GTF -> FASTA
genePredToGtf file nonCano.sorf.formatted.gpe nonCano.sorf.formatted.gtf
bash "$translate_script" nonCano.sorf.formatted.gtf "$fa" $PWD
mv -f prot.fa nonCano.sorf.fa
seqkit fx2tab nonCano.sorf.fa > nonCano.sorf.tab
[[ -s nonCano.sorf.fa && -s nonCano.sorf.tab ]] || die "RibORF sORF fasta/tab not generated."

log "All done. Outputs in: $work_path/$sub_dir/{PRICE,RiboCode,RibORF}"
