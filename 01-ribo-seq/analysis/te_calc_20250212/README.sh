# 计算所有ORF的翻译效率（ORF的RPF reads以及转录本的表达量）

## 1、计算ORF的RPF reads
### 1.1、测试
mkdir -p test && cd test
source activate ribocode
gtf=/home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/sorfs.all.gtf
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output/alignment/SRR2163083_Aligned.sortedByCoord.out.bam
# ORFcount -g <RiboCode_ORFs_result.gtf> -r <ribo-seq genomic mapping file> -f 15 -l 5 -e 100 -m 26 -M 34 -o <ORF.counts>
# not work
ORFcount -g $gtf -r $bam -f 15 -l 5 -e 100 -m 26 -M 34 -o ORF.counts
less $gtf|head -n8|sed 's/gene_id/orf_id/' > test.gtf
# weird，只需要输入ORF所在的gtf，而不需要输入整个转录本的gtf
ORFcount -g test.gtf -r $bam -f 15 -l 5 -e 100 -m 26 -M 34 -o ORF.counts

less $gtf|head -n8 > test.1.gtf
featureCounts -t CDS -g gene_id -a test.1.gtf -o output.counts.txt $bam

p_site_sam=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2015_Science_ChoJ_hippocampus/SRR2163083/output/Ribo_ORFs_add_assemble_20250125/RibORF/corrected.SRR2163083.sam
featureCounts -t CDS -g gene_id -a test.1.gtf -o output.counts.txt.1 $p_site_sam

featureCounts -s 1 -t CDS -g gene_id -a test.1.gtf -o output.counts.txt.2 $p_site_sam
featureCounts -s 2 -t CDS -g gene_id -a test.1.gtf -o output.counts.txt.3 $p_site_sam	
featureCounts -s 2 -t CDS -g gene_id -f -O -a test.1.gtf -o output.counts.txt.4 $p_site_sam

# 只能生成整个文库的质量评估
bam_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/2020_NAR_WangH_E15.5_P42/SRR11218268/output/Ribo-ORFs/PRICE/SRR11218268_Aligned.sortedByCoord.out.bam
samtools index $bam_1
ribotish quality -b $bam_1 -g /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf

### 1.2、对所有的sORF计算RPF reads数目或者p-site的数目
nohup bash Run.cal.p_site.20250213.sh &> log/Run.cal.p_site.20250213.log &
nohup bash Uni.cal.p_site.20250213.sh /home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/essen_info_to_gen_genepred.sorfs.gtf output/ms_sorfs_ribo_reads &> log/Uni.cal.p_site.20250213.log &

## 2、计算ORF的RNA-seq reads

### 2.1、先进行mapping

nohup bash Run.mapping.mouse.rna-seq.sh &> log/Run.mapping.mouse.rna-seq.log &

### 2.2、检查文库类型
# 检查是否为链特异性文库
## 结果表明ERR和SRR216是stranded,SRR5262869和SRR5262868是unstranded,SRR11218257和SRR11218256是reversely stranded
fastq_lst=$PWD/output/qc_rna_seq/filtered_fastq.lst
bam_path=output/mapping_rna_seq_filtered/
output_path=stranded
script=/home/user/BGM/lit/anaconda3/envs/py2/bin/infer_experiment.py
anno_bed=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.bed
find $bam_path -name "*Aligned.sortedByCoord.out.bam" | xargs -I {} sh -c "echo {}; $script -i {} -q 255 -r $anno_bed" > $output_path/whetherStranded.txt

### 2.3、使用featureCount计算reads数目
nohup bash Run.cal.rna_reads.20250213.sh &> log/Run.cal.rna_reads.20250213.log &

### 梳理的时候又发现mapping的时候有一个比对率接近于0，因此前期质控是有必要的
#### 先看下所有fastq的质量
##### 这里脚本写错了，因此直接是把trim以及trim以后的fastqc都跑上了
fastq_lst=output/mapping_rna_seq/rna_fastq.lst
less $fastq_lst| parallel -j 5 --joblog log/fastqc.20250213.log 'bash Uni.qc.mouse.rna_seq.sh {} none output/qc_rna_seq/ yes no no'

cd output/qc_rna_seq
mkdir -p multiqc_tmp_1
find ./ -type f -path "*fastqc/*[0-9]_fastqc*"|xargs -i cp {} ./multiqc_tmp_1
multiqc ./multiqc_tmp_1 --outdir ./ --filename before

mkdir -p multiqc_tmp_2
find ./ -type f -path "*fastqc/*trimmed_fastqc*"|xargs -i cp {} ./multiqc_tmp_2
multiqc ./multiqc_tmp_2 --outdir ./ --filename after

### 2.4、重新mapping
find $PWD/output/qc_rna_seq/ -name "*trimmed.fq.gz" > output/qc_rna_seq/filtered_fastq.lst
nohup bash Uni.mapping.mouse.rna-seq.sh \
	output/mapping_rna_seq_filtered \
	/home/user/data3/lit/resource/genome/mouse/mm39/index/mm39_STARindex \
	$PWD/output/qc_rna_seq/filtered_fastq.lst &> log/Uni.mapping.mouse.rna-seq.log &
#### 计算重新mapping之前以及之后的比对率
find output/mapping_rna_seq_filtered -name "*.final.out" |xargs grep -H "Uniquely mapped reads %" > output/mapping_rna_seq_filtered/Uniquely_mapped_rate.txt
find output/mapping_rna_seq -name "*.final.out" |xargs grep -H "Uniquely mapped reads %" > output/mapping_rna_seq/Uniquely_mapped_rate.txt
# 去掉之前的结果
rm -rf output/mapping_rna_seq
gtf=/home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/sorfs.all.gtf
# 整理
mv output/sorfs_all output/ribo_sorfs
nohup bash Uni.rna_reads.20250219.sh $gtf output/ribo_sorfs/rna_counts &> log/Uni.rna_reads.20250219.ribo_sorfs.log &

# fn的含义是false negative
mv output/ms_sorfs_ribo_reads output/ribo_sorfs_fn
mkdir output/ribo_sorfs_fn/ribo_counts
find output/ribo_sorfs_fn/ -type f |xargs -i cp {} output/ribo_sorfs_fn/ribo_counts
mkdir output/ribo_sorfs/ribo_counts
find output/ribo_sorfs/ -maxdepth 1 -type f |xargs -i mv {} output/ribo_sorfs_fn/ribo_counts

gtf=/home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/sorfs.all.gtf
nohup bash Uni.rna_reads.20250219.sh $gtf output/ribo_sorfs/rna_counts &> log/Uni.rna_reads.20250219.ribo_sorfs.log &
nohup bash Uni.cal.p_site.20250213.sh $gtf output/ribo_sorfs/ribo_counts &> log/Uni.cal.p_site.20250220.ribo_sorfs.log &

gtf_1=/home/user/data3/lit/project/sORFs/03-Cross-anna/gen_genepred/essen_info_to_gen_genepred.sorfs.gtf
nohup bash Uni.rna_reads.20250219.sh $gtf_1 output/ribo_sorfs_fn/rna_counts &> log/Uni.rna_reads.20250219.ribo_sorfs_fn.log &
nohup bash Uni.cal.p_site.20250213.sh $gtf_1 output/ribo_sorfs_fn/ribo_counts &> log/Uni.cal.p_site.20250220.ribo_sorfs_fn.log &

gtf_2=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
nohup bash Uni.rna_reads.20250219.sh $gtf_2 output/anno_orfs/rna_counts &> log/Uni.rna_reads.20250219.anno_orfs.log &
nohup bash Uni.cal.p_site.20250213.sh $gtf_2 output/anno_orfs/ribo_counts &> log/Uni.cal.p_site.20250220.anno_orfs.log &

output_path=$PWD/output/ribo_sorfs_fn
Rscript Uni.Cal_expr_te.R "$output_path/ribo_counts/p_site.counts.txt" "$output_path/rna_counts/" "$output_path/te"

# 得到bam list
find /home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/ -type f -path '*/output/Ribo_ORFs_add_assemble_20250125/RibORF/corrected*sam' > bam_lst/ribo_bam_p_site.lst
find $PWD/output/mapping_rna_seq_filtered/ -name "*Aligned.sortedByCoord.out.bam" > bam_lst/rna_bam_p_site.lst

# 计算libsize，新写了一个uni的脚本
cp /home/user/data2/lit/project/ZNF271/GPB-revision/20240730-calculate-added-data-multi-tissue/calcu-libsize.sh ./Uni.calcu-libsize.sh
## 运行
bash Uni.calcu-libsize.sh bam_lst/rna_bam_p_site.lst libsize/rna
bash Uni.calcu-libsize.sh bam_lst/ribo_bam_p_site.lst libsize/ribo

output_path=$PWD/output/anno_orfs
Rscript Uni.Cal_expr_te.R "$output_path/ribo_counts/p_site.counts.txt" "$output_path/rna_counts/" "./libsize/rna/libsize.txt" "$output_path/te"
Rscript Uni.PCA_analysis.R "./output/anno_orfs/te/rna_seq_rpkm.txt" "./output/anno_orfs/rna_based_pca_corr_rpkm" "rna"
Rscript Uni.PCA_analysis.R "./output/anno_orfs/te/rna_seq_tpm.txt" "./output/anno_orfs/rna_based_pca_corr_tpm" "rna"
Rscript Uni.PCA_analysis.R "./output/anno_orfs/te/TE.txt" "./output/anno_orfs/ribo_based_pca_corr" "ribo"

for output_path in $PWD/output/ribo_sorfs_fn $PWD/output/anno_orfs;do
Rscript Uni.Cal_expr_te.R "$output_path/ribo_counts/p_site.counts.txt" "$output_path/rna_counts/" "./libsize/rna/libsize.txt" "./libsize/ribo/libsize.txt" "$output_path/te"
done

featureCounts -t CDS -g transcript_id -a /home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf -o test/rna_counts_test.txt $bam_lst_1