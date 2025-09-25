#!/usr/bin/sh

################################################
#File Name: README.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Thu 10 Apr 2025 10:26:05 PM CST
################################################

set -eo pipefail

total_cds_bed_1_based=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/sep.add_anno.cds.1-based.bed
total_pep_fa=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/prot.fa
mkdir -p ORFs_bed peptide_fa
for id in $(cat /home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/sep_to_be_checked.txt);do
grep $id $total_cds_bed_1_based > ORFs_bed/$id.ORF.bed
grep -A1 $id $total_pep_fa > peptide_fa/$id.ORF_pep.fa
done

paste <(ls $PWD/ORFs_bed/*bed | sort) <(ls $PWD/peptide_fa/*fa | sort) > input_list.txt
mkdir log
nohup cat input_list.txt | parallel -j 10 --colsep '\t' --joblog log/denovo.20250410.parallel.log bash Uni.check.denovo.20250410.v1.sh {1} {2} &> log/denovo.20250410.log &

path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410
orf_bed_file=$path/ORFs_bed/ENST00000651715.1+chr8:73880420-73880468.ORF.bed 
maf_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/maf
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
python3 $scriptDir/1_extract_multiple_alignments.py \
    -b $orf_bed_file -m $maf_dir -o results_120 -f yes

# 改写脚本，尝试第二步运行
## 如果文件夹下有多个SEP的相关文件，将会一并运行
cd /home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/case
echo "Calculating ancestral sequences and estimnate intact ORF ancestrally"
script=$scriptDir/2.1_ancestral_sequences.v1.sh
script_1=$scriptDir/2.2_parsing_ancestors.v1.py
nwk_file=/home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/tree/120mammal.nwk
bash $script $PWD/ \
    $nwk_file \
    $script_1

script_2=$scriptDir/3_sequence_specificity.v1.py
echo "Performing similarity searches across orthologous regions to trace protein age"
python3 $script_2 $PWD/orfs/ENST00000514422.1+chr4__3589823-3589865.fa $PWD/peptide_fa/ENST00000514422.1+chr4:3589823-3589865.ORF_pep.fa \
    $PWD/orfs/ENST00000514422.1+chr4__3589823-3589865.maf

##### 批量运行步骤二和步骤三 #####
nohup bash Run.step_2_3.sh &> log/step_2_3.log & 

##### 整理步骤二和步骤三的结果 #####
cd results_120/
ls ./prank/*ancestors|egrep -v "nucl|prot"|xargs cat > ancestors.tmp
ls ./orfs/*spec.out|xargs cat > spec.out.tmp
organize(){
grep -v orf_id $1 > tmp.txt
cat <(head -n 1 $1) tmp.txt
rm -rf $1
}
organize ancestors.tmp > ancestors
organize spec.out.tmp > spec.out