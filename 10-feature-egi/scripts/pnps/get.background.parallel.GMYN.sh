#!/bin/bash
set -euo pipefail

#####################
# 基本参数
#####################
# 线程数，默认 30，可以通过环境变量 THREADS 覆盖
THREADS=${THREADS:-30}

# 设置路径
cds_fa_tmp=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.fa
seqkit rmdup -n $cds_fa_tmp > /home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.rmdup.fa
cds_fa=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.rmdup.fa
# head -n 200 /home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.fa > ./test.fa
# cds_fa=$PWD/test.fa
prot_fa=../processed/seve_list_compare/input_for_pnps/all_orfs.cds.noM.fa
out_dir=../processed/pnps/orfs_fa
script=/home/user/data3/lit/project/sORFs/10-feature-egi/scripts/kaks/kaks.bin.sh

mkdir -p "$out_dir"

# 翻译 CDS → 蛋白（如果你已经翻译过，也可以手动注释掉这一行）
faTrans -stop "$cds_fa" "$prot_fa"

# 合并结果文件 + 表头
result_file=$out_dir/all_kaks_results.txt
echo -e "orf_id\tSpecies_pair\tKa\tKs\tKa/Ks\tS-Substitutions\tN-Substitutions\tS-Sites\tN-Sites" > "$result_file"

#####################
# 单个 ORF 的处理函数
#####################
run_kaks_for_orf() {
    local orf_id="$1"

    echo "Processing: $orf_id"

    # 生成安全的文件名
    local safe_id
    safe_id=$(echo "$orf_id" | sed 's/[:|+]/_/g')

    # 提取蛋白序列
    local prot_single="$out_dir/${safe_id}.fa"
    seqkit grep -p "$orf_id" "$prot_fa" > "$prot_single" || true

    if [ ! -s "$prot_single" ]; then
        echo "Warning: No protein sequence found for $orf_id, skipping..."
        return 0
    fi

    # 构建“假 ortholog” 的蛋白 fa
    local prot_ortho="$out_dir/${safe_id}.ortho.prot.fa"
    local prot_ortho_rename="$out_dir/${safe_id}.ortho.prot.rename.fa"
    cat "$prot_single" "$prot_single" > "$prot_ortho"
    seqkit replace -p '.+' -r 'human_{nr}' "$prot_ortho" > "$prot_ortho_rename"

    # 提取 CDS 序列
    local cds_single="$out_dir/${safe_id}.cds.fa"
    seqkit grep -p "$orf_id" "$cds_fa" > "$cds_single"

    # 构建“假 ortholog” 的 CDS fa
    local cds_ortho="$out_dir/${safe_id}.ortho.cds.fa"
    local cds_ortho_rename="$out_dir/${safe_id}.ortho.cds.rename.fa"
    cat "$cds_single" "$cds_single" > "$cds_ortho"
    seqkit replace -p '.+' -r 'human_{nr}' "$cds_ortho" > "$cds_ortho_rename"

    # 子目录 + 跑 KaKs
    local orf_out_dir="$out_dir/$safe_id"
    mkdir -p "$orf_out_dir"

    bash "$script" \
        -p "$prot_ortho_rename" \
        -c "$cds_ortho_rename" \
        -o "$orf_out_dir" \
        -m all \
        -k GMYN

    # 追加结果到合并文件（行顺序可能不再严格按 ORF 保证，但不会缺行）
    if [ -f "$orf_out_dir/kaks.res" ]; then
        awk -v id="$orf_id" 'NR>1{print id"\t"$0}' "$orf_out_dir/kaks.res" >> "$result_file"
    else
        echo "Warning: No kaks.res found for $orf_id"
    fi

    # 如果想清理中间文件，可以取消下面注释
    # rm -f "$prot_single" "$prot_ortho" "$prot_ortho_rename" \
    #       "$cds_single" "$cds_ortho" "$cds_ortho_rename"
}

export -f run_kaks_for_orf
export cds_fa prot_fa out_dir script result_file

#####################
# 并行执行
#####################
# 使用 GNU parallel 并行
seqkit seq -n "$cds_fa" | parallel -j "$THREADS" run_kaks_for_orf {}

# 如果你没有 parallel，也可以用 xargs（把上面 parallel 那行注释掉，换成下面这行）：
# seqkit seq -n "$cds_fa" | xargs -n1 -P "$THREADS" -I{} bash -c 'run_kaks_for_orf "$@"' _ {}

echo "All done! Results saved to: $result_file"
