# 20250822 更改输出路径
output_path=$PWD/02-output-20250821
mkdir -p $output_path
# 20250822 指定gtf路径，而不是从define_annotation_gencode_v41_human中获取
gtf=/home/user/data3/lit/project/sORFs/09-CustomDb/Test_20250801/processed/mkAnno_for_moPepGen/custom.gtf
# bam lst
bam_lst=$PWD/02-output-20250621/featureCounts/bam.lst
##### 3. featureCounts #####
# 首先检查链特异性
mkdir -p $output_path/featureCounts
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o $output_path/featureCounts/rna-seq-counts.txt $(cat $bam_lst)

#### 4.计算表达量 #####
mkdir -p $output_path/expr/
ls $output_path/mapping/*.final.out | xargs -I {} basename {}|sed 's/.R1_Log.final.out//' > $output_path/expr/sample_name.txt
grep 'Uniquely mapped reads number' $output_path/mapping/*.final.out| \
awk '{print $NF}' >  $output_path/expr/lib.txt
paste $output_path/expr/sample_name.txt $output_path/expr/lib.txt | awk '{print $1,$2}' > $output_path/expr/sample_name.lib.txt