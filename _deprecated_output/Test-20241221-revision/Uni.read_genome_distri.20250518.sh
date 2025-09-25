script=/home/user/BGM/lit/anaconda3_bak/envs/py2/bin/read_distribution.py

bam_lst=$1
output_dir=$2

# 默认的 bed 文件路径
default_bed="/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed"
# 检查是否提供了 bed 参数
if [ -z "$3" ]; then
    bed="$default_bed"
else
    bed="$3"
fi
# 打印 bed 文件路径
echo "Using bed file: $bed"

for bam in $(cat $bam_lst);do
sample=$(echo $bam|awk -F'alignment/|_Aligned' '{print $2}')
echo $sample
$script -i $bam -r $bed >  $output_dir/$sample.r_distri.txt
done

cd  $output_dir
out_file="summary_tag_count.tsv"
[ -f $out_file ] && rm -rf $out_file
echo -e "Sample\tType\tTag_count" > $out_file
for file in *r_distri.txt; do
    sample_name=$(basename "$file" .r_distri.txt)
    awk -v sample="$sample_name" 'NR>=6 && NR<=15 {print sample"\t"$1"\t"$3}' "$file" >> $out_file
done