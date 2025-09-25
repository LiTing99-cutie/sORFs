bam_lst=$1
output_path=$2
mkdir -p $output_path
organize_libType_script=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/5.0.Uni.organize.libType.sh
source /home/user/data2/lit/bin/lit_utils.sh
cd $output_path
infer_experiment_human $bam_lst whetherStranded.txt
bash $organize_libType_script whetherStranded.txt
# 得到结果enhanced_results_with_path.tsv