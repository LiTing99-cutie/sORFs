################################################
#File Name: run.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Thu Oct 31 15:08:42 2024
################################################

#!/bin/sh 
for orfs in canonical_orfs denovo_orfs lnc_orfs;do
cut -f 1,10 $workDir/$orfs.txt|tail -n +2|seqkit tab2fx > $workDir/$orfs.fa
done
source activate PrimeNovo-1
workDir=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare
outputDir=$workDir/ISD_per_aa
fastaDir=$workDir/ISD_per_aa
script=/home/user/data3/rbase/opt/aiupred/aiupred.py
aminoAcidFastaDenovo=$workDir/denovo_orfs.fa
aminoAcidFastaPC=$workDir/canonical_orfs.fa
aminoAcidFastaCncRNA=$workDir/lnc_orfs.fa
calculate_average_script=/home/user/data3/rbase/denovo_tumor/denovo_genes/ISD/calculate_average.sh
## 100 de novo genes
aa_types=(aa aa_no_C)
for type in ${aa_types[@]};
do
    echo "#### $type ####"
    echo "--- De novo genes ---"
    sed 's/*$//g' $aminoAcidFastaDenovo | sed 's/C//g' > $fastaDir/100_denovo_genes.$type.fasta
    [ -f $outputDir/100_denovo_genes.$type.ISD.out ] || python $script -i $fastaDir/100_denovo_genes.$type.fasta \
        -o $outputDir/100_denovo_genes.$type.ISD.out -v -g 1
    bash $calculate_average_script $outputDir/100_denovo_genes.$type.ISD.out $workDir/ISD_avg/100_denovo_genes.$type.ISD_avg.txt

    ## PC genes
    echo "--- Protein-coding genes ---"
    sed 's/*$//g' $aminoAcidFastaPC | sed 's/C//g' > $fastaDir/pc_ORFs.shuffle.$type.fa
    [ -f $outputDir/pc_ORFs.shuffle.$type.ISD.out ] || python $script -i $fastaDir/pc_ORFs.shuffle.$type.fa \
        -o $outputDir/pc_ORFs.shuffle.$type.ISD.out -v -g 1
    bash $calculate_average_script $outputDir/pc_ORFs.shuffle.$type.ISD.out $workDir/ISD_avg/pc_ORFs.shuffle.$type.ISD_avg.txt

    ## lncRNA
    echo "--- Cong non-coding RNA genes ---"
    sed 's/*$//g' $aminoAcidFastaCncRNA | sed 's/C//g' > $fastaDir/lncRNA_ORFs.gt30.shuffle.$type.fa
    [ -f $outputDir/lncRNA_ORFs.gt30.shuffle.$type.ISD.out ] || python $script -i $fastaDir/lncRNA_ORFs.gt30.shuffle.$type.fa \
        -o $outputDir/lncRNA_ORFs.gt30.shuffle.$type.ISD.out -v -g 1
    bash $calculate_average_script $outputDir/lncRNA_ORFs.gt30.shuffle.$type.ISD.out $workDir/ISD_avg/lncRNA_ORFs.gt30.shuffle.$type.ISD_avg.txt

done