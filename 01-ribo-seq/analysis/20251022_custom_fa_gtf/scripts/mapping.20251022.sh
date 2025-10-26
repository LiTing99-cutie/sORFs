output_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/scripts/Uni.mapping.using.star.v3.20250725.sh
for clean_fq in $(cat ../processed/clean_fastq.lst);do
bash $script $clean_fq $output_path/alignment
done