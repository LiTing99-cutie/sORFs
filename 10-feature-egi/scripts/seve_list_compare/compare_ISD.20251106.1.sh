#!/bin/bash
################################################
#File Name: run.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Thu Oct 31 15:08:42 2024
################################################

# set -e
# set -u
# set -o pipefail

# 变量定义（移到最前面）
workDir=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare
outputDir=$workDir/ISD_per_aa
fastaDir=$workDir/ISD_per_aa/fasta
script=/home/user/data3/rbase/opt/aiupred/aiupred.py
aminoAcidFastaDenovo=$workDir/denovo_orfs.fa
aminoAcidFastaPC=$workDir/canonical_orfs.fa
aminoAcidFastaCncRNA=$workDir/lnc_orfs.fa
calculate_average_script=/home/user/data3/rbase/denovo_tumor/denovo_genes/ISD/calculate_average.sh

# 创建必要目录
mkdir -p $outputDir $fastaDir $workDir/ISD_avg

# 生成FASTA文件
for orfs in canonical_orfs denovo_orfs lnc_orfs; do
    cut -f 1,10 $workDir/$orfs.txt | tail -n +2 | seqkit tab2fx > $workDir/$orfs.fa
done

# 激活conda环境
source activate PrimeNovo-1

# 主循环
aa_types=(aa aa_no_C)
for type in ${aa_types[@]}; do
    echo "#### $type ####"
    
    # 根据类型决定是否删除C
    if [ "$type" == "aa_no_C" ]; then
        sed_cmd="sed 's/*$//g' | sed 's/C//g'"
    else
        sed_cmd="sed 's/*$//g'"
    fi
    
    ## De novo genes
    echo "--- De novo genes ---"
    eval "$sed_cmd" < $aminoAcidFastaDenovo > $fastaDir/100_denovo_genes.$type.fasta
    if [ ! -f $outputDir/100_denovo_genes.$type.ISD.out ]; then
        python $script -i $fastaDir/100_denovo_genes.$type.fasta \
            -o $outputDir/100_denovo_genes.$type.ISD.out -v -g 1
    fi
    bash $calculate_average_script $outputDir/100_denovo_genes.$type.ISD.out \
        $workDir/ISD_avg/100_denovo_genes.$type.ISD_avg.txt
    
    ## PC genes
    echo "--- Protein-coding genes ---"
    eval "$sed_cmd" < $aminoAcidFastaPC > $fastaDir/pc_ORFs.shuffle.$type.fa
    if [ ! -f $outputDir/pc_ORFs.shuffle.$type.ISD.out ]; then
        python $script -i $fastaDir/pc_ORFs.shuffle.$type.fa \
            -o $outputDir/pc_ORFs.shuffle.$type.ISD.out -v -g 1
    fi
    bash $calculate_average_script $outputDir/pc_ORFs.shuffle.$type.ISD.out \
        $workDir/ISD_avg/pc_ORFs.shuffle.$type.ISD_avg.txt
    
    ## lncRNA
    echo "--- Cong non-coding RNA genes ---"
    eval "$sed_cmd" < $aminoAcidFastaCncRNA > $fastaDir/lncRNA_ORFs.gt30.shuffle.$type.fa
    if [ ! -f $outputDir/lncRNA_ORFs.gt30.shuffle.$type.ISD.out ]; then
        python $script -i $fastaDir/lncRNA_ORFs.gt30.shuffle.$type.fa \
            -o $outputDir/lncRNA_ORFs.gt30.shuffle.$type.ISD.out -v -g 1
    fi
    bash $calculate_average_script $outputDir/lncRNA_ORFs.gt30.shuffle.$type.ISD.out \
        $workDir/ISD_avg/lncRNA_ORFs.gt30.shuffle.$type.ISD_avg.txt
done

echo "All done!"