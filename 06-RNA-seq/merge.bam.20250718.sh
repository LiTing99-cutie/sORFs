merged_bam_dir=/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-20250621/mapping/merge/
mkdir -p $merged_bam_dir
echo "$(date '+%Y-%m-%d %H:%M:%S') Merging Cytoplasm bams"
samtools merge -@ 30 -o $merged_bam_dir/Cytoplasm.bam 02-output-20250621/mapping/*C_*bam
echo "$(date '+%Y-%m-%d %H:%M:%S') Merging Nucleus bams"
samtools merge -@ 30 -o $merged_bam_dir/Nucleus.bam 02-output-20250621/mapping/*N_*bam
echo "$(date '+%Y-%m-%d %H:%M:%S') Merging Total bams"
samtools merge -@ 30 -o $merged_bam_dir/Total.bam $merged_bam_dir/Cytoplasm.bam $merged_bam_dir/Nucleus.bam