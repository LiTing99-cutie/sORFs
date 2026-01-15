source activate biotools

# 输入denovo鉴定的输出路径；按照原来的脚本运行第三步，得到de novo ORF的ortholog
## [以下参数需要调整]
denovo_out_dir=/home/user/data3/rbase/small_peptide/denovo_eval
denovo_prefix=/home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list
bed=$denovo_prefix.bed
orfs_pep_fa=$denovo_prefix.fa
out_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251128_ortholog_generation/results

scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
script_3=$scriptDir/3_sequence_specificity.v3.py
# nucl_dir=$denovo_out_dir/results/orfs/nucl
# pep_dir=$denovo_out_dir/results/orfs/prot
# cp -r $nucl_dir ./
# cp -r $pep_dir ./
nucl_dir=$PWD/nucl
pep_dir=$PWD/prot

# script=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251113_lncORF_denovo_check/scripts/prepare.orf.files.py
# cp $script ./
[ -d ../results/orf_pep/ ] || rm -rf ../results/orf_pep/*
python prepare.orf.files.v2.py \
    "$bed" \
    "$orfs_pep_fa" \
    "$nucl_dir" \
    "$pep_dir" \
    "$out_dir"

# cat $out_dir/input_list.txt |head -n1|parallel -j 1 --colsep '\t' python3 $script_3 {1} {2} {3}
cat $out_dir/input_list.txt | parallel -j 30 --colsep '\t' --joblog ../log/step_3.parallel.log python3 $script_3 {1} {2} {3} &> ../log/step_3.all.log