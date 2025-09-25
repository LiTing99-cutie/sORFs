output_path=$PWD/02-output-isoseq-gtf-20250909
mkdir -p $output_path
gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
# bam lst
bam_lst=$PWD/02-output-20250621/featureCounts/bam.lst
##### 3. featureCounts #####
mkdir -p $output_path/featureCounts
# featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o $output_path/featureCounts/rna-seq-counts.txt $(cat $bam_lst|head -n1)
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o $output_path/featureCounts/rna-seq-counts.txt $(cat $bam_lst)

#### 4.计算表达量 #####
mkdir -p $output_path/expr/
# ls $output_path/mapping/*.final.out | xargs -I {} basename {}|sed 's/.R1_Log.final.out//' > $output_path/expr/sample_name.txt
# grep 'Uniquely mapped reads number' $output_path/mapping/*.final.out| \
# awk '{print $NF}' >  $output_path/expr/lib.txt
# paste $output_path/expr/sample_name.txt $output_path/expr/lib.txt | awk '{print $1,$2}' > $output_path/expr/sample_name.lib.txt