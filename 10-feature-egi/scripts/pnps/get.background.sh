#!/bin/bash

# 设置路径
cds_fa=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.fa
# head -n 200 /home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/all_orfs.cds.fa > ./test.fa
# cds_fa=$PWD/test.fa
prot_fa=../processed/seve_list_compare/input_for_pnps/all_orfs.cds.noM.fa
out_dir=../processed/pnps/orfs_fa
script=/home/user/data3/lit/project/sORFs/10-feature-egi/scripts/kaks/kaks.bin.sh

# 创建输出目录
mkdir -p $out_dir

# 翻译CDS序列（如果还没做的话）
faTrans -stop $cds_fa $prot_fa

# 创建合并结果文件，添加表头
result_file=$out_dir/all_kaks_results.txt
echo -e "orf_id\tSpecies_pair\tKa\tKs\tKa/Ks\tS-Substitutions\tN-Substitutions\tS-Sites\tN-Sites" > $result_file

# 获取所有ORF ID并循环处理
seqkit seq -n $cds_fa | while read orf_id; do
    echo "Processing: $orf_id"
    
    # 创建安全的文件名（替换特殊字符）
    safe_id=$(echo "$orf_id" | sed 's/[:|+]/_/g')
    
    # 提取蛋白序列
    seqkit grep -p "$orf_id" $prot_fa > $out_dir/"$safe_id".fa
    
    # 检查是否成功提取
    if [ ! -s $out_dir/"$safe_id".fa ]; then
        echo "Warning: No protein sequence found for $orf_id, skipping..."
        continue
    fi
    
    # 复制并重命名蛋白序列
    cat $out_dir/"$safe_id".fa $out_dir/"$safe_id".fa > $out_dir/"$safe_id".ortho.prot.fa
    seqkit replace -p '.+' -r 'human_{nr}' $out_dir/"$safe_id".ortho.prot.fa > $out_dir/"$safe_id".ortho.prot.rename.fa
    
    # 提取CDS序列
    seqkit grep -p "$orf_id" $cds_fa > $out_dir/"$safe_id".cds.fa
    
    # 复制并重命名CDS序列
    cat $out_dir/"$safe_id".cds.fa $out_dir/"$safe_id".cds.fa > $out_dir/"$safe_id".ortho.cds.fa
    seqkit replace -p '.+' -r 'human_{nr}' $out_dir/"$safe_id".ortho.cds.fa > $out_dir/"$safe_id".ortho.cds.rename.fa
    
    # 创建输出子目录并运行KaKs分析
    mkdir -p $out_dir/"$safe_id"
    bash $script -p $out_dir/"$safe_id".ortho.prot.rename.fa \
                 -c $out_dir/"$safe_id".ortho.cds.rename.fa \
                 -o $out_dir/"$safe_id" \
                 -m all \
                 -k NG
    
    # 提取结果并追加到合并文件
    if [ -f $out_dir/"$safe_id"/kaks.res ]; then
        awk -v id="$orf_id" 'NR>1{print id"\t"$0}' $out_dir/"$safe_id"/kaks.res >> $result_file
    else
        echo "Warning: No kaks.res found for $orf_id"
    fi
    
    # 可选：清理中间文件以节省空间
    # rm -f $out_dir/"$safe_id".fa $out_dir/"$safe_id".ortho.*.fa
    
done

echo "All done! Results saved to: $result_file"