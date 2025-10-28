database_path=$1
project_path=/rd1/user/lit/project/sORFs
script=$project_path/Uni.Fragpipe.sh
file_hla_1=$project_path/raw_data/HLA/MHC-I/XC03541MHC_P21_HLA1_Slot1-1_1_1100.d
file_hla_2=$project_path/raw_data/HLA/MHC-II/Replicate_2/XC04813MHC_P21_HLA2_Slot1-1_1_18437.d
processed_dir=$(realpath ../processed)
for file in $file_hla_1 $file_hla_2;do
    if [ $file == $file_hla_1 ];then
        output_path=$processed_dir/HLA-I
        workflow=$project_path/fragpipe_config/HLA/Nonspecific-HLA-diaPASEF.workflow
    else
        output_path=$processed_dir/HLA-II
        workflow=$project_path/fragpipe_config/HLA/Nonspecific-HLA-II-diaPASE.workflow
    fi
    mkdir -p $output_path && cd $output_path
    mkdir -p config
    sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > config/fragpipe.workflow
    echo -e "$file\t$(basename $file)\t\tDIA" > config/fragpipe-files.fp-manifest
    bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./ 40
done


