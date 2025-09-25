# 所有的SEP都运行
## peptide_fa的上一级路径；包括了.ORF_pep.fa
work_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/
# work_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/case
## orfs的上一级路径；包括了.maf以及.fa
output_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/results_120
# output_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/case
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
script=$scriptDir/2.1_ancestral_sequences.v1.sh
script_1=$scriptDir/2.2_parsing_ancestors.v1.py
script_2=$scriptDir/3_sequence_specificity.v1.py
nwk_file=/home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/tree/120mammal.nwk

echo "Calculating ancestral sequences and estimnate intact ORF ancestrally"
bash $script $output_path \
    $nwk_file \
    $script_1

echo "Performing similarity searches across orthologous regions to trace protein age"
# 第一个参数为多物种比对的蛋白质文件，第二个参数为人中的SEP，第三个参数为多物种比对的核苷酸文件
# paste <(find $output_path -path "*/orfs/*fa" |sort) \
#  <(find $work_path/peptide_fa -path "*.fa" |sort) \
#  <(find $output_path -path "*/orfs/*maf" |sort) > input_list.1.txt

[ -f input_list.1.txt ] && rm -rf input_list.1.txt
for maf in $(ls $output_path/orfs/*maf);do
ma_pep=$(echo $maf|sed 's/.maf/.fa/')
# ls $ma_pep
id=$(basename -s .maf $maf |sed 's/__/:/')
sep_pep=$work_path/peptide_fa/$id.ORF_pep.fa
# ls $sep_pep
echo -e "$ma_pep\t$sep_pep\t$maf" >> input_list.1.txt
done
cat input_list.1.txt | parallel -j 10 --colsep '\t' --joblog log/step_2_3.parallel.log python3 $script_2 {1} {2} {3}
