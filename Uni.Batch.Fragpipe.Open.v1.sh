#!/usr/bin/sh

################################################
#File Name: Uni.Batch.Fragpipe.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月14日 星期四 10时54分22秒
################################################

set -eo pipefail
# 20250206 增加可选参数null_set，用于指定null_set；默认为nocleavage；增加多种消化酶
all_expr=$1
metadata=$2
database_path=$3
output_path=$4
# nonspecific or nocleavage
null_set=$5
script=/rd1/user/lit/project/sORFs/Uni.Fragpipe.sh

mkdir -p $output_path && cd $output_path

for i in $(awk -F, '{print $2}' $metadata | tail -n +2 | sort | uniq); do
    echo "***Processing $i..."
    mkdir -p $i && cd $i
    awk -F, -v OFS=',' -v i="$i" '$2 == i' $metadata > MS_metadata_$i.csv
    grep -f <(awk -F, '{print "_" $1 "_"}' MS_metadata_$i.csv) $all_expr > all_expr_$i.txt
    mkdir -p config

    ## 生成manifest文件
    less all_expr_$i.txt | xargs -I {} sh -c 'echo -e "{}\t$(basename {})\t\tDDA"' > config/fragpipe-files.fp-manifest

    ## 生成workflow文件
    workflow_path=/rd1/user/lit/project/sORFs/config/Open
    if [ "$i" == "Null" ]; then
        if [ "$null_set" == "nonspecific" ]; then
            workflow=$workflow_path/Nonspecific_6_50_prot_1.workflow
        else
            workflow=$workflow_path/Nocleavage_6_50_prot_1.workflow
        fi
    elif [ "$i" == "Trypsin" ]; then
        workflow=$workflow_path/Trypsin_6_50_prot_1.workflow
    elif [ "$i" == "Trypsin_LysC" ]; then
        workflow=$workflow_path/Trypsin_LysC_6_50_prot_1.workflow
    elif [ "$i" == "Trypsin_LysN" ]; then
        workflow=$workflow_path/Trypsin_LysN_6_50_prot_1.workflow
    elif [ "$i" == "Trypsin_Chymotrypsin" ]; then
        workflow=$workflow_path/Trypsin_Chymotrypsin_6_50_prot_1.workflow
    elif [ "$i" == "Trypsin_ArgC" ]; then
        workflow=$workflow_path/Trypsin_ArgC_6_50_prot_1.workflow
    elif [ "$i" == "Trypsin_AspN" ]; then
        workflow=$workflow_path/Trypsin_AspN_6_50_prot_1.workflow
    elif [ "$i" == "Trypsin_GluC" ]; then
        workflow=$workflow_path/Trypsin_GluC_6_50_prot_1.workflow
    else
        echo "Warning: No specific workflow found for $i. Skipping."
        continue
    fi

    ## 替换数据库路径并生成配置文件
    sed 's|database.db-path=.*|database.db-path='"$database_path"'|' $workflow > config/fragpipe.workflow

    ## 运行fragpipe
    time bash $script config/fragpipe.workflow config/fragpipe-files.fp-manifest ./ 40

    ## 返回上级目录
    cd ..
done