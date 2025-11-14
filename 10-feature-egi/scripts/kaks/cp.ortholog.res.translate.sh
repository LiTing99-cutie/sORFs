#!/bin/bash

id_file="/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/canonical_orfs.lnc_orfs.intergenic_orfs.id.txt"
denovo_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251113_lncORF_denovo_check
target_dir="../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_cano_lnc_intergenic_sample"

mkdir -p "$target_dir"

echo "Step 1: 并行拷贝nucl文件..."

export denovo_dir target_dir

copy_nucl() {
    orig_id="$1"
    # safe_id=$(echo "$orig_id" | sed 's/[^A-Za-z0-9._-]/__/g')
    safe_id=$(echo "$orig_id" | sed 's/[^A-Za-z0-9._-]/__/g' | sed 's/__\+/__/g')
    
    nucl_file=$(find "$denovo_dir" -type f -name "${safe_id}.*ortholog.nucl.fa" | head -1)
    
    if [ -s "$nucl_file" ]; then
        cp "$nucl_file" "$target_dir/"
        echo "✓ $(basename "$nucl_file")"
    elif [ -f "$nucl_file" ]; then
        echo "✗ 文件为空: $(basename "$nucl_file") (0字节)" >&2
    else
        echo "✗ 文件不存在: $orig_id" >&2
    fi
}

export -f copy_nucl

cat "$id_file" |head -n50| parallel -j 50 copy_nucl
# cat "$id_file" | parallel -j 50 copy_nucl

echo ""
echo "Step 2: 批量翻译nucl → pep..."

# 批量翻译所有nucl文件
translate_one() {
    nucl_file="$1"
    pep_file="${nucl_file%.nucl.fa}.pep.fa"
    
    if transeq -sequence "$nucl_file" -outseq "$pep_file" -frame 1 -table 1 -clean 2>/dev/null; then
        echo "✓ 翻译: $(basename "$pep_file")"
    else
        echo "✗ 翻译失败: $(basename "$nucl_file")" >&2
    fi
}

export -f translate_one

# find "$target_dir" -name "*.ortholog.nucl.fa" | head -n50|parallel -j 50 translate_one
find "$target_dir" -name "*.ortholog.nucl.fa" |parallel -j 50 translate_one

echo ""
echo "Done!"
echo "目标目录: $target_dir"
echo "nucl文件数: $(ls -1 "$target_dir"/*.nucl.fa 2>/dev/null | wc -l)"
echo "pep文件数: $(ls -1 "$target_dir"/*.pep.fa 2>/dev/null | wc -l)"