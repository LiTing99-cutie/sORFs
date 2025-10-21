# chr=chr1
chr=$1
# chr=chrY
proj_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251013_denovo_check/processed
## peptide_fa的上一级路径；包括了.ORF_pep.fa
work_path=$proj_path/get_input/per_chr/peptide_fa/$chr
## orfs的上一级路径；包括了.maf以及.fa
output_path=$proj_path/check_denovo/$chr
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
# script_2_1=$scriptDir/2.1_ancestral_sequences.v1.sh
script_2_1=$scriptDir/2.1_ancestral_sequences.v2.sh
script_2_2=$scriptDir/2.2_parsing_ancestors.v1.py
# script_3=$scriptDir/3_sequence_specificity.v1.py
script_3=$scriptDir/3_sequence_specificity.v2.py
nwk_file=/home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/tree/120mammal.nwk

##### 第二步 ##### 
echo "Calculating ancestral sequences and estimnate intact ORF ancestrally"
mkdir -p ../log/$chr
bash $script_2_1 -j 80 $output_path \
    $nwk_file \
    $script_2_2 &> ../log/$chr/run.step_2.log

echo "Performing similarity searches across orthologous regions to trace protein age"
# 第一个参数为多物种比对的蛋白质序列，第二个参数为人中ORF的蛋白质序列，第三个参数为多物种比对的核苷酸序列

##### 第三步【独立于上一步】##### 
map_tsv="../processed/work_mafsInRegion/$chr/$chr.name_map.tsv"
# 生成 SAFE->ORIG 的去重映射
awk -v OFS='\t' '
  {
    split($1, a, "__exon"); safe=a[1];
    if (!(safe in seen)) { print safe, $2; seen[safe]=1 }
  }
' "$map_tsv" > ../processed/work_mafsInRegion/$chr/$chr.orf_safe_map.tsv
map2=../processed/work_mafsInRegion/$chr/$chr.orf_safe_map.tsv

: > input_list.txt
while IFS= read -r -d '' maf; do
  base=$(basename "$maf" .maf)            # SAFE 名
  ma_pep="${maf%.maf}.fa"                 # 同前缀 .fa
  # 从映射表拿到 original id
  orig=$(awk -v s="$base" '$1==s{print $2; exit}' "$map2")
  ORF_pep="$work_path/${orig}.ORF_pep.fa"
  printf "%s\t%s\t%s\n" "$ma_pep" "$ORF_pep" "$maf" >> input_list.txt
done < <(find "$output_path/orfs" -maxdepth 1 -type f -name '*.maf' -print0)

# 需要更改ORF_pep的文件名，再进行输入
cat input_list.txt | parallel -j 30 --colsep '\t' --joblog ../log/$chr/step_3.parallel.log python3 $script_3 {1} {2} {3} &> ../log/$chr/step_3.all.log

##### 整理第二步和第三步的结果 #####
cd $output_path
ls ./prank/*ancestors|egrep -v "nucl|prot"|xargs cat > ancestors.tmp
ls ./orfs/*spec.out|xargs cat > spec.out.tmp
organize(){
grep -v orf_id $1 > tmp.txt
cat <(head -n 1 $1) tmp.txt
rm -rf $1 tmp.txt
}
organize ancestors.tmp > ancestors
organize spec.out.tmp > spec.out