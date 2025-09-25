/home/user/data/lit/project/ZNF271/data/ribo-seq/2020-nature-E-MTAB-7247/human_ribo_seq/fastq/
/home/user/data/lit/project/ZNF271/data/ribo-seq/2020-nature-E-MTAB-7247/human_rna_seq/fastq/

mkdir -p 2020-nature-E-MTAB-7247/ribo 2020-nature-E-MTAB-7247/RNA
ln -s /home/user/data/lit/project/ZNF271/data/ribo-seq/2020-nature-E-MTAB-7247/human_ribo_seq/fastq/* 2020-nature-E-MTAB-7247/ribo
ln -s /home/user/data/lit/project/ZNF271/data/ribo-seq/2020-nature-E-MTAB-7247/human_rna_seq/fastq/* 2020-nature-E-MTAB-7247/RNA

mkdir -p 2022-MC-GSE182377/ribo 2022-MC-GSE182377/RNA
ln -s /home/user/data/lit/project/ZNF271/data/ribo-seq/2022-MC-GSE182377/human/thalamus/* 2022-MC-GSE182377/ribo
ln -s /home/user/data/xieyn/peptidome/data/human/2022_Molecular_Cell/RNA-seq/* 2022-MC-GSE182377/RNA

mkdir -p 2022-NN/ribo 2022-NN/RNA
ln -s /home/user/data/rbase/ribosome_profiling/Duffy_2022_NatNeurosci-brain/fastq/*fastq.gz 2022-NN/ribo/
ln -s /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/2022-NN/fastq/*/*.fastq.gz 2022-NN/RNA
ln -s /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/2022-NN/metadata/SraRunTable.txt 2022-NN/

mkdir -p 2014-JN-GSE51424/ribo 2014-JN-GSE51424/RNA
ln -s /home/user/data/lit/project/ZNF271/data/ribo-seq/2014-JN-GSE51424/human/normal_brain/*fastq.gz 2014-JN-GSE51424/ribo/
ln -s /home/user/data/lit/project/ZNF271/data/ribo-seq/2014-JN-GSE51424/human/normal_brain/rna-seq/*fastq.gz 2014-JN-GSE51424/RNA/
# qc on previous runned rna-seq data
export script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/Uni.qc.mouse.rna_seq.sh
export output_path=/home/user/data3/lit/project/ZNF271/data/rna-seq/brain/2022-NN/qc
mkdir -p $output_path/log
ls /home/user/data3/lit/project/ZNF271/data/rna-seq/brain/2022-NN/fastq/*/*.fastq.gz| \
parallel -j 5 --joblog $output_path/log/qc.20250217.log 'bash $script {} none $output_path no yes yes'
pushd $output_path && mkdir -p multiqc && find ./ -path "*/fastqc/*" -o -name ""|xargs -i cp {} ./multiqc && multiqc ./multiqc && popd


cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/te_calc_20250212/Uni.qc.mouse.rna_seq.sh Uni.qc.single.sh

# fastqc on all ribo-seq files
bash Uni.qc.batch.sh /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc \
<(find $PWD -path "*/ribo/*fastq.gz") \
yes no no
## 有些下载的质量不佳，需要重新下载，重新run
bash Uni.qc.batch.sh /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc \
<(find $PWD -path "*/ribo/human_brain_ribo_1.fastq.gz") \
yes no no

outdir=$PWD/qc
pushd /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc 
[ -d $outdir/multiqc_before ] || mkdir -p $outdir/multiqc_before
find ./ -path "*/fastqc/*" |xargs -i cp {} $outdir/multiqc_before && multiqc --outdir $outdir -n multiqc_before.html $outdir/multiqc_before && popd

# fastqc on remaining rna-seq files
bash Uni.qc.batch.sh /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc_rna_seq \
<(find $PWD -path "*/2020*/RNA/*fastq.gz" -o -path "*/2022-MC*/RNA/*fastq.gz") \
yes no no

## 有些下载的质量不佳，需要重新下载，重新run
bash Uni.qc.batch.sh /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc_rna_seq \
<(find $PWD -path "*/2020*/RNA/human_brain_rna_2.fastq.gz" -o -path "*/2022-MC*/RNA/SRR15513227_2.fastq.gz") \
yes no no

outdir=$PWD/qc_rna_seq
pushd /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc_rna_seq 
[ -d $outdir/multiqc_before ] || mkdir -p $outdir/multiqc_before
find ./ -path "*/fastqc/*" |xargs -i cp {} $outdir/multiqc_before && multiqc --outdir $outdir -n multiqc_before.html $outdir/multiqc_before && popd

# 针对adaptor质量不合格的RNA-seq进行adaptor trimming [remaining]
bash Uni.qc.batch.sh /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/qc_rna_seq \
<(find $PWD -path "*/2020*/RNA/*fastq.gz") \
no yes yes

##### 测试 ##### 
fq=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/2022-NN/ribo/SRR15906460.fastq.gz
# universal+8个A
## 去不干净
bash Uni.qc.single.mul.ada.sh $fq \
	$PWD/adaptor.fa \
	$PWD/tmp \
	no yes yes
mv tmp/SRR15906460/ tmp/SRR15906460_1

# 只有8个A
## 可以去干净
bash Uni.qc.single.mul.ada.sh $fq \
	$PWD/adaptor.1.fa \
	$PWD/tmp \
	no yes yes

# universal+--polyA选项
## 报错
bash Uni.qc.single.mul.ada.polyA.sh $fq \
	$PWD/adaptor.2.fa \
	$PWD/tmp_2 \
	no yes yes

# cutadapt+10个A
cutadapt -a AAAAAAAAAA -q 20,20 -m 20 -u 3 -u -2 \
         --cores 30 -o tmp/tmp.filtered.fastq.gz \
         $fq > cutadapt.filter.log
fastqc tmp/tmp.filtered.fastq.gz -O tmp

# 10个A+5端去除3碱基(work)
bash Uni.qc.single.2022_NN.sh $fq \
	AAAAAAAAAA \
	$PWD/tmp_3 \
	no yes yes
##### 测试 END ##### 

# 针对adaptor质量不合格的Ribo-seq进行adaptor trimming（分研究进行）
## 2022-NN
bash Uni.qc.batch.v1.sh $PWD/qc/2022_NN \
	<(find $PWD -path "*/2022-NN/ribo/*fastq.gz") \
	no yes yes \
	Uni.qc.single.2022_NN.sh \
	AAAAAAAAAA

## 2022-MC
bash Uni.qc.batch.v1.sh $PWD/qc/2022-MC \
	<(find $PWD -path "*/2022-MC*/ribo/*fastq.gz") \
	no yes yes \
	Uni.qc.single.sh \
	AGATCGGAAGAGCACACGTCT
### 测试一下别的adaptor
fq=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Human_organized/2022-MC-GSE182377/ribo/SRR15513149_GSM5527674_Brain_2_RiboSeq_Homo_sapiens_RNA-Seq.fastq.gz
adapter=/home/user/data/rbase/ribosome_profiling/Chothani_2022_MolCell-brain_embryonic/fastq/adapters.fa
bash Uni.qc.single.mul.ada.sh $fq \
	$adapter \
	$PWD/tmp_3 \
	no yes yes

## 2020-Nature
bash Uni.qc.batch.v1.sh $PWD/qc/2020-Nature \
	<(find $PWD -path "*/2020-nature*/ribo/*fastq.gz") \
	no yes yes \
	Uni.qc.single.sh \
	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

## 2014
bash Uni.qc.batch.v1.sh $PWD/qc/2014-JN-GSE51424 \
	<(find $PWD -path "*/2014-JN-GSE51424/ribo/*fastq.gz") \
	yes yes yes \
	Uni.qc.single.sh \
	CTGTAGGCACCATCAAT

##### 新添加的RNA-seq数据跑一下fastqc
bash Uni.qc.batch.sh $PWD/qc_rna_seq/2014-JN-GSE51424 \
<(find $PWD -path "*/2014-JN-GSE51424/RNA/*fastq.gz") \
yes no no

### 20250227 重新排列文件
mkdir qc_rna_seq/2020-Nature
mv qc_rna_seq/human_brain* qc_rna_seq/2020-Nature