bam_lst=$1
output_path=$2
mkdir -p $output_path
organize_libType_script=$PWD/Uni.organize.libType.sh
source /home/user/data2/lit/bin/lit_utils.sh
cd $output_path
infer_experiment_human $bam_lst whetherStranded.txt
bash $organize_libType_script whetherStranded.txt
# 得到结果enhanced_results_with_path.tsv