out_dir=../processed/evolution_para/kaks/denovo
mkdir -p $out_dir
denovo_dir_1=$(realpath ../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_denovo)
# 随机选取的cano,lnc以及intergenic【需要起始密码子翻译成对应密码子】
bash ./kaks/kaks.parallel.sh $denovo_dir_1 $out_dir &> ../log/kaks.parallel.denovo.log

# 合并，增加多余信息
merge_ksks(){
    base_dir=$1

    # 提取表头
    find $base_dir -name "*aln.kaks" | head -1 | xargs head -1 > header.txt
    awk '{print "ORF_ID\t"$0}' header.txt > $base_dir/all_kaks_merged.txt

    # 并行处理
    find $base_dir -name "*aln.kaks" | \
    parallel -j 8 '
        orf_id=$(echo {} | grep -oP "PB\.[^/]+")
        tail -n +2 {} | awk -v orf="$orf_id" "{print orf\"\t\"\$0}"
    ' >> $base_dir/all_kaks_merged.txt

    echo "完成！"
}

merge_ksks $out_dir
