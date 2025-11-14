source activate biotools

SCRIPT_DIR=/home/user/data3/rbase/small_peptide/denovo_eval/scripts
maf_dir="/home/user/data3/lit/project/sORFs/05-denovo-status/maf"
out_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/results
tmp_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/tmp
bed=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/canonical.intergenic_orfs.lnc_orfs.bed
nohup python $SCRIPT_DIR/extract_multiple_alignment.py \
--input_bed $bed \
--maf_dir $maf_dir \
--output_dir $out_dir \
--work_dir $tmp_dir \
--focal hg38 \
--start_codon ATG,CTG,GTG,TTG,ACG \
--cov_thresh 0.95 \
--max_procs 24 &

##### 第三步【独立于上一步】##### 
# # 需要更改ORF_pep的文件名，再进行输入
# # "$ma_pep" "$ORF_pep" "$maf"
nucl_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/results/orfs/nucl
pep_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/results/orfs/prot
orf_pep_dir=$(realpath ../results/orf_pep)
orfs_pep_fa=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/canonical.intergenic_orfs.lnc_orfs.fa
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
script_3=$scriptDir/3_sequence_specificity.v3.py
# # 生成ORF_pep.fa
# mkdir -p $orf_pep_dir
# cut -f 4 $bed|sort -u > $orf_pep_dir/orf.id.txt
# # 第一步：生成原始ID到safe_id的映射
# awk '{
#     original = $0
#     safe = $0
#     gsub(/[^A-Za-z0-9._-]+/, "__", safe)
#     print original "\t" safe
# }' $orf_pep_dir/orf.id.txt > $orf_pep_dir/id_mapping.txt
# # 第二步：使用映射文件进行grep
# python generate_orf_pep_fa.py ${orf_pep_dir} $orf_pep_dir/id_mapping.txt $orfs_pep_fa

# [ -f $out_dir/input_list.txt ] && rm -rf $out_dir/input_list.txt
# find ${nucl_dir} -name "*.nucl.fa" | while read fa; do
#     id=$(basename "$fa" .nucl.fa)
#     prot=$(ls ${pep_dir}/chr*/${id}.prot.fa)
#     orf=$(ls ${orf_pep_dir}/${id}.ORF_pep.fa)
#     nucl=$(ls ${nucl_dir}/chr*/${id}.nucl.fa)
#     echo -e "$prot\t$orf\t$nucl"
# done >> $out_dir/input_list.txt
bed=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/canonical.intergenic_orfs.lnc_orfs.bed
orfs_pep_fa=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/canonical.intergenic_orfs.lnc_orfs.fa
nucl_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/results/orfs/nucl
pep_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/results/orfs/prot
out_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/2021113_lncORF_denovo_check/results

python prepare.orf.files.py \
    "$bed" \
    "$orfs_pep_fa" \
    "$nucl_dir" \
    "$pep_dir" \
    "$out_dir"

# cat $out_dir/input_list.txt |head -n1|parallel -j 1 --colsep '\t' python3 $script_3 {1} {2} {3}
nohup cat $out_dir/input_list.txt | parallel -j 30 --colsep '\t' --joblog ../log/step_3.parallel.log python3 $script_3 {1} {2} {3} &> ../log/step_3.all.log &