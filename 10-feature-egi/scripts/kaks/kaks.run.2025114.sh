out_dir=../processed/evolution_para/kaks
mkdir -p $out_dir/cano_lnc_intergenic_noncano_sampled
mkdir -p $out_dir/noncano_sampled
denovo_dir_1=$(realpath ../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_cano_lnc_intergenic_sample)
# 随机选取的cano,lnc以及intergenic【需要起始密码子翻译成对应密码子】
bash ./kaks/kaks.parallel.sh $denovo_dir_1 $out_dir/cano_lnc_intergenic_noncano_sampled &> ../log/kaks.parallel.1.log
# 随机选取的noncano
denovo_dir_2=$(realpath ../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_noncano_sample)
# bash ./kaks/kaks.parallel.sh $denovo_dir_2 $out_dir/noncano_sampled &> ../log/kaks.parallel.2.log

# 合并，增加多余信息
merge_ksks(){
    base_dir=$1

    # 提取表头
    find $base_dir -name "*aln.kaks" | head -1 | xargs head -1 > header.txt
    awk '{print "ORF_ID\t"$0}' header.txt > $base_dir/all_kaks_merged.txt

    # 并行处理
    find $base_dir -name "*aln.kaks" | \
    parallel -j 8 '
        orf_id=$(basename $(dirname $(dirname {})))
        tail -n +2 {} | awk -v orf="$orf_id" "{print orf\"\t\"\$0}"
    ' >> $base_dir/all_kaks_merged.txt

    echo "完成！"
}
# orf_id=$(echo {} | grep -oP "PB\.[^/]+")
merge_ksks /home/user/data3/lit/project/sORFs/10-feature-egi/processed/evolution_para/kaks/cano_lnc_intergenic_noncano_sampled
merge_ksks /home/user/data3/lit/project/sORFs/10-feature-egi/processed/evolution_para/kaks/noncano_sampled