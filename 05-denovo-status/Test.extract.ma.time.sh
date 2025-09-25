source activate denovo
work_dir=$PWD
maf_dir=/home/user/data/rbase/120_mammal_alignment/maf
orf_bed_file=$work_dir/ORFs_bed/ENSG00000100433.ORF.bed
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
cd test_time/plan_1
mkdir -p tmp
# 304m1.601s 5h
time python3 $scriptDir/1_extract_multiple_alignments.py \
    -b $orf_bed_file -m $maf_dir -o results_120 -f yes 

cd ../..
cp /home/user/data/rbase/120_mammal_alignment/maf/chr14.maf ./
mkdir -p test_time/plan_2 && cd test_time/plan_2
mkdir -p tmp
# 14 min
time python3 $scriptDir/1_extract_multiple_alignments.py \
    -b $orf_bed_file -m /home/user/data3/lit/project/sORFs/05-denovo-status -o results_120 -f yes 

# 移动maf到data3目录下
nohup find /home/user/data/rbase/120_mammal_alignment/maf/ -type f -name "chr*.maf" | grep -E '/chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\.maf$'|xargs -I {} cp {} ./maf &