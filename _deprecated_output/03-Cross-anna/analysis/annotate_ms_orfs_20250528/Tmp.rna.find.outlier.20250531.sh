
mapping_res_dir=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble
output_dir=$PWD/output/S2/rna_seq_qc
mkdir -p $output_dir
##### 0. 得到样本的不同输出list #####
find $mapping_res_dir -name *_Log.final.out > $output_dir/Log.final.out.lst
##### 2. 唯一比对数量以及唯一比对率 #####
source /home/user/data3/lit/project/sORFs/sORFs.utils.sh
extract_info $output_dir/Log.final.out.lst "Uniquely mapped reads %" \
    "$output_dir/Uniquely_mapped_reads_rate.txt"
extract_info $output_dir/Log.final.out.lst "Uniquely mapped reads number" \
    "$output_dir/Uniquely_mapped_reads_number.txt"
merge_files_by_sample "$output_dir/Uniquely_mapped_reads_rate.txt" "$output_dir/Uniquely_mapped_reads_number.txt" \
	"$output_dir/Uniquely_mapped_reads_rate_number.txt"