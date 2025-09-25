# 注释掉的包使用conda install -c bioconda * 安装
conda create --name denovo --file ./Denovo_genes-tumors/evolution_orfs/requirements.txt
mkdir -p ORFs_bed peptide_fa
cp /home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/ORFs_bed/ENSG00000100433.ORF.bed ./ORFs_bed
cp /home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/peptide_fa/ENSG00000100433.ORF_pep.fa ./peptide_fa
work_dir=$PWD
maf_dir=/home/user/data/rbase/120_mammal_alignment/maf
nwk_file=/home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/tree/120mammal.nwk
orf_bed_file=$work_dir/ORFs_bed/ENSG00000100433.ORF.bed
pep_fasta_file=$work_dir/peptide_fa/ENSG00000100433.ORF_pep.fa
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
mkdir -p $work_dir/tmp/results_120
mkdir -p $work_dir/results_120
        
echo "Calculating multiple alignments. Recommend to parallelize the searches for multiple ORFs."
# python3 $scriptDir/1_extract_multiple_alignments.py \
#     -b $orf_bed_file -m $maf_dir -o results_120 -f yes

echo "Calculating ancestral sequences and estimate intact ORF ancestrally"
bash $scriptDir/2.1_ancestral_sequences.sh $work_dir/results_120 \
    $nwk_file \
    $orf_bed_file \
    $pep_fasta_file \
	$scriptDir/2.2_parsing_ancestors.py

echo "Performing similarity searches across orthologous regions to trace protein age"
python3 $scriptDir/3_sequence_specificity.py --prot_dir $work_dir/results_120 --prot_tar $pep_fasta_file

echo "Curating these de novo ORFs manually using primate synteny data"
### maltiple sequence alignments for these ORFs and define "common disablers"