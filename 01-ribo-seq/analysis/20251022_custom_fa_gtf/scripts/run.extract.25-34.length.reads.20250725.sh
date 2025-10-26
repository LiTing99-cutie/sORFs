bam_dir=$1
output_dir=$2
mkdir -p $output_dir
find $bam_dir -name "*sort*bam" > $output_dir/bam.lst
for bam in $(cat $output_dir/bam.lst);do
    echo "Processing $bam"
    bash Uni.extract_specific_length_reads.v1.20250520.sh $bam $output_dir/$(basename $bam)
done