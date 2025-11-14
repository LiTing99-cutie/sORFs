out_dir=../processed/evolution_para/kaks
mkdir -p $out_dir/cano_lnc_intergenic_noncano_sampled
mkdir -p $out_dir/noncano_sampled
denovo_dir_1=$(realpath ../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_cano_lnc_intergenic_sample)
# 随机选取的cano,lnc以及intergenic【需要起始密码子翻译成对应密码子】
bash ./kaks/kaks.parallel.sh $denovo_dir_1 $out_dir/cano_lnc_intergenic_noncano_sampled &> ../log/kaks.parallel.1.log
# 随机选取的cano
denovo_dir_2=$(realpath ../processed/seve_list_compare/input_for_ortholog_extraction/ortholog_noncano_sample)
bash ./kaks/kaks.parallel.sh $denovo_dir_2 $out_dir/noncano_sampled &> ../log/kaks.parallel.2.log
