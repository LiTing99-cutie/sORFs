id_file=../processed/seve_list_compare/r1_m_minus_non_cano_orfs.sample.2.5k.lt.150aa.id.txt
denovo_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251029_denovo_check/processed
target_dir="../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_noncano_r1_m_minus_sample"

mkdir -p "$target_dir"

# 导出函数和变量供parallel使用
export denovo_dir target_dir

sanitize_and_copy() {
    orig_id="$1"
    # safe_id=$(echo "$orig_id" | sed 's/[^A-Za-z0-9._-]/__/g')
    safe_id=$(echo "$orig_id" | sed 's/[^A-Za-z0-9._-]/__/g' | sed 's/__\+/__/g')
    
    pep_file=$(find "$denovo_dir" -type f -name "${safe_id}.*ortholog.pep.fa" | head -1)
    nucl_file=$(find "$denovo_dir" -type f -name "${safe_id}.*ortholog.nucl.fa" | head -1)
    
    [ -n "$pep_file" ] && cp "$pep_file" "$target_dir/" && echo "✓ $(basename "$pep_file")"
    [ -n "$nucl_file" ] && cp "$nucl_file" "$target_dir/" && echo "✓ $(basename "$nucl_file")"
}

export -f sanitize_and_copy

# 并行处理 (8个任务同时运行)
cat "$id_file" | parallel -j 50 sanitize_and_copy

echo "Done!"