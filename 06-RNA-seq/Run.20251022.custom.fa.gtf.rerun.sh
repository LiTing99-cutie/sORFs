output_path=$PWD/02-output-custom-fa-gtf-20251022
mkdir -p $output_path
gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
rawdata_path=/home/user/data3/lit/project/sORFs/06-RNA-seq/01-rawdata
##### 2. mapping #####
# mkdir -p $output_path/mapping
# mkdir -p log
# # 从第二个开始运行
# find $rawdata_path -name "*fastq.gz"|grep -v p21_C_3 > $output_path/mapping/fastq.lst
# # 更换脚本和STAR运行的设置
# bash bin/1.2.Uni.rna.mapping.assemble.human.20251023.sh\
# 	$output_path/mapping/fastq.lst \
# 	$output_path/mapping/ &> log/mapping.$(date '+%Y%m%d').log

##### 3. featureCounts #####
mkdir -p $output_path/featureCounts
ls $output_path/mapping/*_Aligned.sortedByCoord.out.bam > $output_path/featureCounts/bam.lst
bam_lst=$output_path/featureCounts/bam.lst
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o $output_path/featureCounts/rna-seq-counts.txt $(cat $bam_lst)

#### 4.计算表达量 #####
mkdir -p $output_path/expr/
ls $output_path/mapping/*.final.out | xargs -I {} basename {}|sed 's/.R1_Log.final.out//' > $output_path/expr/sample_name.txt
grep 'Uniquely mapped reads number' $output_path/mapping/*.final.out| \
awk '{print $NF}' >  $output_path/expr/lib.txt
paste $output_path/expr/sample_name.txt $output_path/expr/lib.txt | awk '{print $1,$2}' > $output_path/expr/sample_name.lib.txt