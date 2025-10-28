# 此脚本用于计算mzML格式的质谱下机数据的ms2谱图数量（其中.d格式的质谱数据通过fragpipe转换为mzML格式）
# 输入的mzML list
mzML_lst=$1
# 输出的ms2_count文件
output_file=$2
# 清空输出文件
: > $PWD/tmp.1.txt
: > $PWD/tmp.2.txt
cat $mzML_lst | xargs -I {} basename {} | grep -oP '(?<=min_)(.*?)(?=_Slot)' >> $PWD/tmp.1.txt
cat $mzML_lst | xargs -I {} grep -c "spectrum index=" {} >> $PWD/tmp.2.txt
paste $PWD/tmp.1.txt $PWD/tmp.2.txt > $output_file
rm -rf $PWD/tmp.1.txt $PWD/tmp.2.txt