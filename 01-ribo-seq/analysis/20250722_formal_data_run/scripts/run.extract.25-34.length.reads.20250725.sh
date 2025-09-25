bam_dir=$(realpath "../processed/batch_1/alignment_star_20250725_change_index_para/")
output_dir=$(realpath "../processed/batch_1/filtered_bam")
mkdir -p $output_dir
find $bam_dir -name "*sort*bam" > $output_dir/bam.lst
for bam in $(cat $output_dir/bam.lst);do
    echo "Processing $bam"
    bash Uni.extract_specific_length_reads.v1.20250520.sh $bam $output_dir/$(basename $bam)
done