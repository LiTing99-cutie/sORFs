#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# -------------------- 默认参数 --------------------
CONDA_ENV="base"
NPROC=8
BY_SAMPLE_DIR="../processed/by_sample"
LOGROOT="/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/log"
RESULTS_DIR="../results"

RPF_MIN="0.2"
PS_MIN="0.1"
MIN_C="0.2"

PEP_ORF_MERGED="../MS_res_from_Galaxy/pep.orf.merged.txt"
INTENSITY_MERGED="../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv"
file2=../MS_res_from_Galaxy/merged.peptide.tsv                            # 列: Peptide_I_L_equal Source Source_number Sample

RUN_EACH_SAMPLE="./run.all.steps.20250918.v3.sh"   # <- 新版每样本脚本名
SPLIT_BY_SAMPLE="split_by_sample.py"
MERGE_ADD_SAMPLE="merge_add_sample.py"
ORF_FOLD_AND_OCC="orf_fold_and_occurrence.py"
ORF_MEAN_REL_IBAQ="orf_mean_relative_ibaq.py"
MERGE_TABLES_ALL="merge_tables_all_sample.py"
SPEARMAN_V2="spearman_corr.v2.py"

RPF="/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt"
ISO="/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt"
RNA="/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt"
GENE_ANNO="/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt"

ALL_IBAQ_B=""   # 运行时根据 RESULTS_DIR 生成
ALL_ASSIGN=""
FOLDED=""
UNIQ_PEP=""
MEAN_REL_IBAQ=""
MERGED_FINAL=""
SPEARMAN_OUT=""

usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  --by-sample-dir DIR     默认: $BY_SAMPLE_DIR
  --logroot DIR           默认: $LOGROOT
  --results-dir DIR       默认: $RESULTS_DIR
  --nproc N               默认: $NPROC
  --conda-env NAME        默认: $CONDA_ENV

  --rpf-min FLOAT         默认: $RPF_MIN
  --ps-min  FLOAT         默认: $PS_MIN
  --min-c   FLOAT         默认: $MIN_C

示例：
  $0 --by-sample-dir ../processed/by_sample --logroot ./log \\
     --results-dir ../results --nproc 8 --rpf-min 0.2 --ps-min 0.1 --min-c 0.2
EOF
  exit 1
}

# -------------------- 解析参数 --------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --by-sample-dir) BY_SAMPLE_DIR="$2"; shift 2;;
    --logroot)       LOGROOT="$2"; shift 2;;
    --results-dir)   RESULTS_DIR="$2"; shift 2;;
    --nproc)         NPROC="$2"; shift 2;;
    --conda-env)     CONDA_ENV="$2"; shift 2;;
    --rpf-min)       RPF_MIN="$2"; shift 2;;
    --ps-min)        PS_MIN="$2"; shift 2;;
    --min-c)         MIN_C="$2"; shift 2;;
    -h|--help)       usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

log(){ printf '[%(%F %T)T] %s\n' -1 "$*"; }
ensure(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 不在 PATH 中"; exit 127; }; }
check_file(){ [[ -s "$1" ]] || { echo "ERROR: 缺失或空文件：$1"; exit 2; }; }

main() {
  mkdir -p "$BY_SAMPLE_DIR" "$RESULTS_DIR" "$LOGROOT"

  ALL_IBAQ_B="${RESULTS_DIR}/all_samples.ibaq_b_with_total.tsv"
  ALL_ASSIGN="${RESULTS_DIR}/all_samples.assignments.unique.post.tsv"
  FOLDED="${RESULTS_DIR}/ibaq_orf_folded.tsv"
  UNIQ_PEP="${RESULTS_DIR}/unique_peptide_counts.tsv"
  UNIQ_PEP_with_source="${RESULTS_DIR}/unique_peptide_counts_with_source.tsv"
  MEAN_REL_IBAQ="${RESULTS_DIR}/ibaq_orf_mean_relative.tsv"
  MERGED_FINAL="${RESULTS_DIR}/orfs_merged_final.tsv"
  SPEARMAN_OUT="${RESULTS_DIR}/spearman_results.tsv"

  # 0) conda
  log "Activate conda env: $CONDA_ENV"
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"

  # 依赖
  for t in python3 awk sort xargs; do ensure "$t"; done
  for py in "$SPLIT_BY_SAMPLE" "$MERGE_ADD_SAMPLE" "$ORF_FOLD_AND_OCC" "$ORF_MEAN_REL_IBAQ" "$MERGE_TABLES_ALL" "$SPEARMAN_V2"; do
    [[ -f "$py" ]] || { echo "ERROR: 脚本不存在：$py"; exit 3; }
  done
  check_file "$PEP_ORF_MERGED"
  check_file "$INTENSITY_MERGED"

  # 1) 按样本拆分
  SAMPLE_LIST="$BY_SAMPLE_DIR/sample_list.txt"
  log "Split by sample -> $BY_SAMPLE_DIR"
  python3 "$SPLIT_BY_SAMPLE" \
    --pep_orf_merged "$PEP_ORF_MERGED" \
    --intensity_merged "$INTENSITY_MERGED" \
    --out_base "$BY_SAMPLE_DIR" \
    --sample_list "$SAMPLE_LIST"
  check_file "$SAMPLE_LIST"

  # 2) 并行每样本（无 nohup，阻塞直到完成）
  log "Run per-sample (P=$NPROC) with thresholds: rpf_min=$RPF_MIN ps_min=$PS_MIN min_c=$MIN_C"
  export RUN_EACH_SAMPLE BY_SAMPLE_DIR LOGROOT RPF_MIN PS_MIN MIN_C
  cat "$SAMPLE_LIST" \
  | xargs -P "$NPROC" -I{} bash -lc '
      sample="{}"
      d="'"$LOGROOT"'/'"{}"'/log"
      mkdir -p "$d"
      bash "'"$RUN_EACH_SAMPLE"'" --sample "$sample" \
           --by-sample-dir "'"$BY_SAMPLE_DIR"'" \
           --rpf-min "'"$RPF_MIN"'" --ps-min "'"$PS_MIN"'" --min-c "'"$MIN_C"'" \
        > "$d/driver.1.log" 2>&1
    '

  # 3) 合并 iBAQ
  log "Merge iBAQ -> $ALL_IBAQ_B"
  python3 "$MERGE_ADD_SAMPLE" \
    --root "$BY_SAMPLE_DIR" \
    --name "ibaq_b_with_total.tsv" \
    --out "$ALL_IBAQ_B" \
    --mode dirnameN \
    --depth 2

  # 4) 折叠 ORF
  log "Fold ORFs -> $FOLDED"
  python3 "$ORF_FOLD_AND_OCC" --in "$ALL_IBAQ_B" --out "$FOLDED"

  # 5) 合并 assignments
  log "Merge assignments -> $ALL_ASSIGN"
  python3 "$MERGE_ADD_SAMPLE" \
    --root "$BY_SAMPLE_DIR" \
    --name "assignments.unique.post.tsv" \
    --out "$ALL_ASSIGN" \
    --mode dirnameN \
    --depth 2

  # 6) 统计 unique peptide
  log "Count unique peptides -> $UNIQ_PEP"
  awk -F'\t' '
    NR==1{ if(tolower($1) ~ /^peptide/ || tolower($2) ~ /^assigned/) next }
    { p=$1; id=$2; if(p==""||id=="") next
      key=p "\t" id; if(!(key in seen)){ seen[key]=1; cnt[id]++ } }
    END{ print "ORF_id\tUnique_peptide_n"
         for(id in cnt) print id "\t" cnt[id] }
  ' "$ALL_ASSIGN" | sort -k1,1 > "$UNIQ_PEP"

  file1=$ALL_ASSIGN   # 列: Peptide Assigned_protein_id reason Sample

  awk -F'\t' -v OFS='\t' '
    # ---------- 读 文件2：建立 (Peptide_I_L_equal,Sample) 是否含 msfragger_closed 的映射 ----------
    NR==FNR {
      if (FNR==1) { for(i=1;i<=NF;i++) h2[$i]=i; next }
      k = $h2["Peptide_I_L_equal"] "\t" $h2["Sample"]
      if ($h2["Source"] ~ /msfragger_closed/) closed[k]=1
      next
    }

    # ---------- 读 文件1：统计 ----------
    FNR==1 { for(i=1;i<=NF;i++) h1[$i]=i; next }

    {
      p = $h1["Peptide"]
      id= $h1["Assigned_protein_id"]   # ORF_id
      s = $h1["Sample"]
      if (p=="" || id=="") next

      pair = p "\t" id
      # 总 unique（跨样本去重）
      if (!(pair in seen)) { seen[pair]=1; cnt[id]++ }

      # 若该 (Peptide, Sample) 在文件2中标注过 msfragger_closed，则记为 closed-unique
      k = p "\t" s
      if (k in closed && !(pair in seen_closed)) { seen_closed[pair]=1; cnt_closed[id]++ }
    }

    END {
      print "ORF_id","Unique_peptide_n","Unique_peptide_n_msfragger_closed"
      for (id in cnt) {
        c2 = (id in cnt_closed) ? cnt_closed[id] : 0
        print id, cnt[id], c2
      }
    }
  ' "$file2" "$file1" | sort -k1,1 > "$UNIQ_PEP_with_source"

  echo "[OK] 写出：$UNIQ_PEP_with_source"

  # 7) 计算相对 iBAQ
  log "Mean relative iBAQ -> $MEAN_REL_IBAQ"
  python3 "$ORF_MEAN_REL_IBAQ" --in "$ALL_IBAQ_B" --metric iBAQ_A --out "$MEAN_REL_IBAQ"

  # 8) 合并全部
  for f in "$RPF" "$ISO" "$RNA" "$GENE_ANNO" "$UNIQ_PEP_with_source" "$MEAN_REL_IBAQ"; do check_file "$f"; done
  log "Merge all -> $MERGED_FINAL"
  python3 "$MERGE_TABLES_ALL" \
    --ibaq "$FOLDED" --rpf "$RPF" --iso "$ISO" --rna "$RNA" \
    --gene_anno "$GENE_ANNO" --uniq_pep "$UNIQ_PEP_with_source" --rel_ibaq "$MEAN_REL_IBAQ" \
    --out "$MERGED_FINAL"

  # 9) 相关性
  log "Spearman -> $SPEARMAN_OUT"
  python3 "$SPEARMAN_V2" \
    --merged "$MERGED_FINAL" \
    --out "$SPEARMAN_OUT" \
    --pairs "A:RPF_RPKM" "A:mean_relative_iBAQ" "RPF_RPKM:mean_relative_iBAQ"

  log "All done."
}

main "$@"
