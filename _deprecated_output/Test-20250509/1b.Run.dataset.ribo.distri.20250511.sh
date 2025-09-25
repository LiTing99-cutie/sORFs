script=/home/user/BGM/lit/anaconda3/envs/py2/bin/read_distribution.py
bed=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.bed
for bam in $(cat $output_dir/Aligned.sortedByCoord.out.bam.lst);do
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