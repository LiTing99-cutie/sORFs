output_path=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S5
rna_seq_bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/merged.sorted.bam

mkdir -p $output_path

infer_experiment_human(){
bam_lst=$1
# whetherStranded.txt
output_file=$2
script=/home/user/BGM/lit/anaconda3/envs/py2/bin/infer_experiment.py
anno_bed=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed
less $bam_lst | xargs -I {} sh -c "echo {}; $script -i {} -q 255 -r $anno_bed -s 2000000"  > $output_file
}

bam_lst=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/bam.lst
infer_experiment_human $bam_lst ./whetherStranded.txt

nohup bash 5.1.merge_bam_count_reads.sh & 
