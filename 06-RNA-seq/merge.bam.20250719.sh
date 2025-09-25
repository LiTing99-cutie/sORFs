merged_bam_dir=/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-20250621/mapping/merge/
[ -f $merged_bam_dir/Total.bam.bai ] || samtools index -@ 30 $merged_bam_dir/Total.bam
[ -f $merged_bam_dir/Cytoplasm.bam.bai ] || samtools index -@ 30 $merged_bam_dir/Cytoplasm.bam
[ -f $merged_bam_dir/Nucleus.bam.bai ] || samtools index -@ 30 $merged_bam_dir/Nucleus.bam
samtools view -@ 30 -c $merged_bam_dir/Cytoplasm.bam > $merged_bam_dir/Cytoplasm.bam.count
samtools view -@ 30 -c $merged_bam_dir/Nucleus.bam > $merged_bam_dir/Nucleus.bam.count
