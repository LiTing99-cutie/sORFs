mkdir bin
mkdir -p 02-output/fastqc
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/Uni*.sh bin
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/1.2.Uni.rna.mapping.assemble.human.20250227.sh bin
##### 1. fastqc #####
# 1 output_dir; 2 bam lst; 3 fastqc; 4 rm adapter; 5 fastqc; 6 script; 7 adapter sequence; all is necessary
bash bin/Uni.qc.batch.v1.sh 02-output/fastqc \
<(find $PWD/01-rawdata/MJ20250407316-MJ-R-20250418027/cleandata/ -type f -name "*fastq.gz"|grep "/L") \
yes no no $PWD/bin/Uni.qc.single.sh "xx"

##### 2. mapping #####
mkdir -p 02-output/mapping
mkdir log
ls $PWD/01-rawdata/MJ20250407316-MJ-R-20250418027/cleandata/*fastq.gz > 02-output/mapping/fastq.lst
nohup bash bin/1.2.Uni.rna.mapping.assemble.human.20250508.sh \
	$PWD/02-output/mapping/fastq.lst \
	$PWD/02-output/mapping/ &> log/mapping.20250508.log &

##### 3. featureCounts #####
# 首先检查链特异性
mkdir 02-output/featureCounts
ls $PWD/02-output/mapping/*_Aligned.sortedByCoord.out.bam > 02-output/featureCounts/bam.lst
bash /home/user/data2/lit/bin/lit_utils.sh
infer_experiment_human 02-output/featureCounts/bam.lst 02-output/featureCounts/whetherStranded.txt 
define_annotation_gencode_v41_human
## featureCounts的链特异性参数-s需要根据特定的文库来修改
featureCounts -s 2 -p --countReadPairs -T 10 -a $gtf -o 02-output/featureCounts/rna-seq-counts.txt $(cat 02-output/featureCounts/bam.lst)

#### 4.计算表达量 #####
mkdir 02-output/expr/
ls /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/mapping/*.final.out | xargs -I {} basename {}|sed 's/.R1_Log.final.out//' > 02-output/expr/sample_name.txt
grep 'Uniquely mapped reads number' /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/mapping/*.final.out| \
awk '{print $NF}' >  02-output/expr/lib.txt
paste 02-output/expr/sample_name.txt 02-output/expr/lib.txt | awk '{print $1,$2/2}' > 02-output/expr/sample_name.lib.txt
