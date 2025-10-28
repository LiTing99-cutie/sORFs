#!/usr/bin/sh

################################################
#File Name: Gene_enzyme.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年12月03日 星期二 16时11分13秒
################################################

set -eo pipefail

# Modify the enzyme settings in a workflow file
# 20241218 add the no digestion enzyme set
# v1 20250114 add the Trypsin_ArgC,Trypsin_AspN,Trypsin_GluC

enzyme=$1
output_file=$2
# 基于哪个workflow进行修改
workflow=$3

if [ "$enzyme" == "Trypsin_LysN" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=lysn
    cut_1=KR
    cut_2=K
    no_cut_1=P
    no_cut_2=
    sense_1=C
    sense_2=N
elif [ "$enzyme" == "Trypsin" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=null
    cut_1=KR
    cut_2=
    no_cut_1=P
    no_cut_2=
    sense_1=C
    sense_2=C
elif [ "$enzyme" == "Trypsin_LysC" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=lysc
    cut_1=KR
    cut_2=K
    no_cut_1=P
    no_cut_2=P
    sense_1=C
    sense_2=C
elif [ "$enzyme" == "Trypsin_ArgC" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=argc
    cut_1=KR
    cut_2=R
    no_cut_1=P
    no_cut_2=P
    sense_1=C
    sense_2=C
elif [ "$enzyme" == "Trypsin_AspN" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=aspn
    cut_1=KR
    cut_2=D
    no_cut_1=P
    no_cut_2=
    sense_1=C
    sense_2=N
elif [ "$enzyme" == "Trypsin_GluC" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=gluc
    cut_1=KR
    cut_2=DE
    no_cut_1=P
    no_cut_2=P
    sense_1=C
    sense_2=C
elif [ "$enzyme" == "Trypsin_Chymotrypsin" ]; then
    enzyme_name_1=trypsin
    enzyme_name_2=chymotrypsin
    cut_1=KR
    cut_2=FLWY
    no_cut_1=P
    no_cut_2=P
    sense_1=C
    sense_2=C
elif [ "$enzyme" == "Nonspecific" ]; then
    enzyme_name_1=nonspecific
    enzyme_name_2=nonspecific
    cut_1=-
    cut_2=-
    no_cut_1=
    no_cut_2=
    sense_1=C
    sense_2=C
elif [ "$enzyme" == "Nocleavage" ]; then
    enzyme_name_1=nocleavage
    enzyme_name_2=nocleavage
    cut_1=-
    cut_2=-
    no_cut_1=
    no_cut_2=
    sense_1=C
    sense_2=C
else
    echo "Warning: No specific workflow found for $enzyme. Skipping."
    continue
fi

sed -e 's|msfragger.misc.fragger.enzyme-dropdown-1=.*|msfragger.misc.fragger.enzyme-dropdown-1='"$enzyme_name_1"'|' \
    -e 's|msfragger.misc.fragger.enzyme-dropdown-2=.*|msfragger.misc.fragger.enzyme-dropdown-2='"$enzyme_name_2"'|' \
    -e 's|msfragger.search_enzyme_name_1=.*|msfragger.search_enzyme_name_1='"$enzyme_name_1"'|' \
    -e 's|msfragger.search_enzyme_name_2=.*|msfragger.search_enzyme_name_2='"$enzyme_name_2"'|' \
    -e 's|msfragger.search_enzyme_cut_1=.*|msfragger.search_enzyme_cut_1='"$cut_1"'|' \
    -e 's|msfragger.search_enzyme_cut_2=.*|msfragger.search_enzyme_cut_2='"$cut_2"'|' \
    -e 's|msfragger.search_enzyme_nocut_1=.*|msfragger.search_enzyme_nocut_1='"$no_cut_1"'|' \
    -e 's|msfragger.search_enzyme_nocut_2=.*|msfragger.search_enzyme_nocut_2='"$no_cut_2"'|' \
    -e 's|msfragger.search_enzyme_sense_1=.*|msfragger.search_enzyme_sense_1='"$sense_1"'|' \
    -e 's|msfragger.search_enzyme_sense_2=.*|msfragger.search_enzyme_sense_2='"$sense_2"'|' \
    $workflow > $output_file

if [ "$enzyme" == "Trypsin_Chymotrypsin" ]; then
    sed -i 's|msfragger.allowed_missed_cleavage_2=.*|msfragger.allowed_missed_cleavage_2=4|' $output_file
fi