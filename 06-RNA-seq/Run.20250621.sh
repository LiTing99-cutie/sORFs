output_path=$PWD/02-output-20250621
mkdir -p $output_path

##### 1. fastqc #####
mkdir -p $output_path/fastqc
# 1 output_dir; 2 bam lst; 3 fastqc; 4 rm adapter; 5 fastqc; 6 script; 7 adapter sequence; all is necessary
bash bin/Uni.qc.batch.v1.sh $output_path/fastqc \
<(find $PWD/01-rawdata/organized_20250624 -name "*fastq.gz") \
yes no no $PWD/bin/Uni.qc.single.sh "xx"

##### 2. mapping #####
mkdir -p $output_path/mapping
mkdir -p log
find $PWD/01-rawdata/organized_20250624 -name "*fastq.gz" > $output_path/mapping/fastq.lst
bash bin/1.2.Uni.rna.mapping.assemble.human.20250508.sh \
	$output_path/mapping/fastq.lst \
	$output_path/mapping/ &> log/mapping.$(date '+%Y%m%d').log

##### 3. featureCounts #####
# 首先检查链特异性
mkdir -p $output_path/featureCounts
ls $output_path/mapping/*_Aligned.sortedByCoord.out.bam > $output_path/featureCounts/bam.lst
source /home/user/data2/lit/bin/lit_utils.sh
infer_experiment_human $output_path/featureCounts/bam.lst $output_path/featureCounts/whetherStranded.txt 
define_annotation_gencode_v41_human
## featureCounts的链特异性参数-s需要根据特定的文库来修改
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o $output_path/featureCounts/rna-seq-counts.txt $(cat $output_path/featureCounts/bam.lst)

#### 4.计算表达量 #####
mkdir -p $output_path/expr/
ls $output_path/mapping/*.final.out | xargs -I {} basename {}|sed 's/.R1_Log.final.out//' > $output_path/expr/sample_name.txt
grep 'Uniquely mapped reads number' $output_path/mapping/*.final.out| \
awk '{print $NF}' >  $output_path/expr/lib.txt
paste $output_path/expr/sample_name.txt $output_path/expr/lib.txt | awk '{print $1,$2}' > $output_path/expr/sample_name.lib.txt
