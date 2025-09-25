#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SAMPLE=""
BY_SAMPLE_DIR=""
RPF_MIN="0.2"
PS_MIN="0.1"
MIN_C="0.2"

usage(){
  cat <<EOF
Usage: $0 --sample NAME --by-sample-dir DIR [--rpf-min F --ps-min F --min-c F]

示例：
  $0 --sample 21pcw_1_3_30K_LC_T --by-sample-dir ../processed/by_sample \\
     --rpf-min 0.2 --ps-min 0.1 --min-c 0.2
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --by-sample-dir) BY_SAMPLE_DIR="$2"; shift 2;;
    --rpf-min) RPF_MIN="$2"; shift 2;;
    --ps-min)  PS_MIN="$2"; shift 2;;
    --min-c)   MIN_C="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -n "$SAMPLE" && -n "$BY_SAMPLE_DIR" ]] || usage

proj_path="/home/user/data3/lit/project/sORFs"
analysis_path="/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map"

id_map="$proj_path/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/id.map.txt"
orf_seq_len="$proj_path/09-CustomDb/formal_20250821/results/custom_db_20250826/orf_seq_len.tsv"

# 本样本 I/O
output_dir="$BY_SAMPLE_DIR/$SAMPLE"
mkdir -p "$output_dir/pick_orf" "$output_dir/pep_assign" "$output_dir/quant"

pep_orf="$output_dir/pep.orf.txt"
peptide_intensity_IL="$output_dir/peptide_intensity_IL.tsv"

echo "[`date '+%F %T'`] Start id.convert.20250911.sh"
bash id.convert.20250911.sh "$pep_orf" "$id_map" \
     "$output_dir/pick_orf/pep.orf.rePicked.txt"

echo "[`date '+%F %T'`] Start pep_to_protein_assign.20250915.alt.col.py"
python pep_to_protein_assign.20250915.alt.col.py \
  --pep2prot "$output_dir/pick_orf/pep.orf.rePicked.txt" \
  --isoform   "$proj_path/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt" \
  --orf_psites "$proj_path/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt" \
  --nca        "$proj_path/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt" \
  --rpf_metric RPF_codon_coverage --rpf_min "$RPF_MIN" \
  --ps_metric  Psites_codon_coverage --ps_min "$PS_MIN" \
  --min_c "$MIN_C" \
  --outdir "$output_dir/pep_assign"

echo "[`date '+%F %T'`] Start make_orf_tables.py"
python make_orf_tables.py \
  --cand_filtered "$output_dir/pep_assign/candidates.filtered.tsv" \
  --orf_seq_len   "$orf_seq_len" \
  --isoform_info  "$proj_path/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt" \
  --out_detail    "$output_dir/pep_assign/orf_detail.noCanonSharing.tsv" \
  --out_orf       "$output_dir/pep_assign/orf_folded_counts.tsv"

echo "[`date '+%F %T'`] Start count_theoretical_peptides.py"
python3 count_theoretical_peptides.py \
  --in "$output_dir/pep_assign/orf_folded_counts.tsv" \
  --sample "$SAMPLE" \
  --out "$output_dir/pep_assign/orf_folded_counts_with_theo_pep.tsv"

echo "[`date '+%F %T'`] Start compute_ibaq_len_norm.v1.py"
python compute_ibaq_len_norm.v1.py \
  --file1 "$output_dir/pep_assign/assignments.unique.post.tsv" \
  --file2 "$peptide_intensity_IL" \
  --file3 "$output_dir/pep_assign/orf_folded_counts_with_theo_pep.tsv" \
  --file3-cols all \
  -o "$output_dir/quant/ibaq_b_with_total.tsv"

echo "[`date '+%F %T'`] Done sample: $SAMPLE"
