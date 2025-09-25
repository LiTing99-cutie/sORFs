################################################
#File Name: run.ribo-seq.analysis.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Fri 03 Nov 2023 10:19:44 AM CST
################################################

#!/bin/sh 

##并发运行脚本，并控制并发数
# 设置并发的进程数
thread_num=5
a=$(date +%H%M%S)
# mkfifo
tempfifo="my_temp_fifo"
mkfifo ${tempfifo}
# 使文件描述符为非阻塞式
exec 6<>${tempfifo}
rm -f ${tempfifo}

# 为文件描述符创建占位信息
for ((i=1;i<=${thread_num};i++))
do
{
    echo 
}
done >&6 #事实上就是在fd6中放置了$thread个回车符

WORK_DIR=/home/user/data3/rbase/denovo_tumor/denovo_genes/ribosome_profiling/reanalyzed/Chothani_2022_MolCell-brain_embryonic
DATA_DIR=/home/user/data/rbase/ribosome_profiling/Chothani_2022_MolCell-brain_embryonic/fastq
BAM_DIR=/home/user/data/rbase/ribosome_profiling/Chothani_2022_MolCell-brain_embryonic/bam
TRIMMOMATIC_DIR=/home/user/data3/rbase/opt/Trimmomatic-0.39
GENOME=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/fasta/Homo_sapiens.GRCh38.primary_assembly.genome
GPE=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/Homo_sapiens.GRCh38.gencode.v43.annotation.gpe
GTF=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/Homo_sapiens.GRCh38.gencode.v43.annotation.gtf
DENOVO_HC_DIR=/home/user/data3/rbase/denovo_tumor/denovo_genes/high-confident
DENOVO_HC_GTF=$DENOVO_HC_DIR/gtf/denovo_genes.high_confident.hominoid.gtf
RIBO=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/fasta/ribosome_rna.cDNA
RiboTISH_DIR=$WORK_DIR/RiboTISH
CandidateORF_DIR=/home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/candidateORF
INDEX=/home/user/data3/rbase/denovo_tumor/normal_rna_seq/GTEx/star_index

samples=(SRR155131{48..58})

candidates=(ENSG00000172927 ENSG00000176912 ENSG00000196542 ENSG00000236081 ENSG00000261175 ENSG00000198547 ENSG00000203930 ENSG00000205704)

# 1. fastqc control quality
echo "### fastqc control quality ###"
for sample in ${samples[@]};
do
    # 一个read -u6命令执行一次，就从FD6中减去一个回车符，然后向下执行
    # 当FD6中没有回车符时，就停止，从而实现线程数量控制
    read -u6
    {
        [ -d $WORK_DIR/fastqc ] || mkdir -p $WORK_DIR/fastqc

        # fastqc qulity control
        echo "-- for $sample --"
        [ -f $WORK_DIR/fastqc/${sample}_fastqc.html ] || fastqc -o $WORK_DIR/fastqc -t 10 $DATA_DIR/${sample}.fastq.gz

        # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
        echo >&6
    }&
done
wait

## merge
[ -f $WORK_DIR/fastqc/multiqc_report.html ] || multiqc $WORK_DIR/fastqc/ --outdir $WORK_DIR/fastqc/

# 2. fastp filter low-quality reads and adapters
echo "### fastp filter low-quality reads and adapters ###"
for sample in ${samples[@]};
do
    read -u6
    {
        [ -d $DATA_DIR/${sample} ] || mkdir -p $DATA_DIR/${sample}
        cd $DATA_DIR/${sample}
        # run fastp
        echo "-- for $sample --"
        # TRALING: remove low quality bases from the end
        # HEADCROP: remove the specified number of bases from the beginning of the read
        # MINLEN: remove the fall below the specific minimal length
        [ -f ${sample}.filtered.fastq.gz ] || \
        java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar SE -phred33 $DATA_DIR/${sample}.fastq.gz ${sample}.filtered.fastq.gz \
            ILLUMINACLIP:../adapters.fa:2:30:5 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20 -threads 30 > trimmomatic.filter.log
        # --trim_poly_x:  enable polyX trimming in 3' ends.
        # --length_required: reads shorter than length_required will be discarded, default is 15. (int [=15]) Default 0 means no limitation.
        # [ -f ${sample}.filtered.fastq.gz ] || \
        # fastp -i $DATA_DIR/${sample}.fastq.gz -o ${sample}.filtered.fastq.gz \
        #     --trim_poly_x --length_required 20 \
        #     -w 10 -j fastp.fliter.json -h fastp.filter.html \
        #     2> fastp.filter.log
        # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
        echo >&6
    }&
done
wait

## fastqc again of clear fastq
echo "### fastqc control quality of clear fastq ###"
for sample in ${samples[@]};
do
    # 一个read -u6命令执行一次，就从FD6中减去一个回车符，然后向下执行
    # 当FD6中没有回车符时，就停止，从而实现线程数量控制
    read -u6
    {
        [ -d $WORK_DIR/fastqc_clear ] || mkdir -p $WORK_DIR/fastqc_clear

        # fastqc qulity control
        echo "-- for $sample --"
        [ -f $WORK_DIR/fastqc_clear/${sample}.filtered_fastqc.html ] || \
        fastqc -o $WORK_DIR/fastqc_clear -t 10 $DATA_DIR/${sample}/${sample}.filtered.fastq.gz
        # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
        echo >&6
    }&
done
wait
## merge
[ -f $WORK_DIR/fastqc_clear/multiqc_report.html ] || multiqc $WORK_DIR/fastqc_clear/ --outdir $WORK_DIR/fastqc_clear/


# 3. remove ribosome RNA
# index 
echo "### index for  bowtie ###"
# bowtie2-build $RIBO.fa $RIBO
# align
echo "### remove ribosome RNA by bowtie ###"
for sample in ${samples[@]};
do
    read -u6
    {
        cd $DATA_DIR/${sample}
        echo "-- for $sample --"
        # -U Files with unpaired reads. Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
        # -S File for SAM output (default: stdout)
        # --un <path>        write unpaired reads that didn't align to <path>
        # --al <path>        write unpaired reads that aligned at least once to <path>
        [ -f ${sample}.filtered.unaligned.fq.gz ] || bowtie2 -q --phred33 -N 1 -x $RIBO \
            -U ${sample}.filtered.fastq.gz -S ${sample}.filtered.ribo.sam --un-gz ${sample}.filtered.unaligned.fq.gz \
            --al-gz ${sample}.filtered.aligned.fq.gz 1> bowtie.report.log 2>&1
        [ -f ${sample}.filtered.ribo.sam ] && rm ${sample}.filtered.ribo.sam
        # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
        echo >&6
    }&
done
wait

# 4. aligning and mapping
echo "### aligning and mapping by STAR ###"
for sample in ${samples[@]};
do
    if [ ! -f $BAM_DIR/${sample}/${sample}.uniq.sorted.bam.flagstate ]; then
        read -u6
        {
            echo "-- for $sample --"
            [ -d $BAM_DIR/${sample} ] || mkdir -p $BAM_DIR/${sample}
            cd $BAM_DIR/${sample}
            ### Alignment 1st Pass.
            time STAR \
            --genomeDir $INDEX \
            --readFilesIn $DATA_DIR/${sample}/${sample}.filtered.unaligned.fq.gz \
            --readFilesCommand zcat \
            --alignEndsType EndToEnd \
            --runThreadN 20 \
            --seedSearchStartLmax 15 \
            --outFilterMismatchNmax 2 \
            --outSJfilterOverhangMin 30 8 8 8 \
            --outFilterScoreMin 0 \
            --outFilterScoreMinOverLread 0.66 \
            --outFilterMatchNmin 0 \
            --outFilterMatchNminOverLread 0.66 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            --outFileNamePrefix ./${sample}.
            
            # unique
            time samtools view -h ${sample}.Aligned.sortedByCoord.out.bam | grep -P '^@|NH:i:1\t' | \
                samtools sort -@ 20 -m 5G -o ${sample}.uniq.sorted.bam
            
            # index
            samtools index -@ 20 ${sample}.Aligned.sortedByCoord.out.bam
            # index
            samtools index -@ 20 ${sample}.uniq.sorted.bam
            
            # stat
            samtools flagstat -@ 20 ${sample}.Aligned.sortedByCoord.out.bam > ${sample}.Aligned.sortedByCoord.out.bam.flagstate
            samtools flagstat -@ 20 ${sample}.uniq.sorted.bam > ${sample}.uniq.sorted.bam.flagstate
            # remove
            [ -f ${sample}.Aligned.sortedByCoord.out.bam ] && rm ${sample}.Aligned.sortedByCoord.out.bam

            # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
            echo >&6
        }&
    fi
done
wait

### merge bam for replicates ###
echo "### merge bam ###"
mkdir -p $BAM_DIR/SRR15513148_49_50_51_52
samtools merge -@ 20 -o $BAM_DIR/SRR15513148_49_50_51_52/SRR15513148_49_50_51_52.uniq.sorted.bam \
    $BAM_DIR/SRR15513148/SRR15513148.uniq.sorted.bam \
    $BAM_DIR/SRR15513149/SRR15513149.uniq.sorted.bam \
    $BAM_DIR/SRR15513150/SRR15513150.uniq.sorted.bam \
    $BAM_DIR/SRR15513151/SRR15513151.uniq.sorted.bam \
    $BAM_DIR/SRR15513152/SRR15513152.uniq.sorted.bam &
mkdir -p $BAM_DIR/SRR15513153_54_55_56_57_58
samtools merge -@ 20 -o $BAM_DIR/SRR15513153_54_55_56_57_58/SRR15513153_54_55_56_57_58.uniq.sorted.bam \
    $BAM_DIR/SRR15513153/SRR15513153.uniq.sorted.bam \
    $BAM_DIR/SRR15513154/SRR15513154.uniq.sorted.bam \
    $BAM_DIR/SRR15513155/SRR15513155.uniq.sorted.bam \
    $BAM_DIR/SRR15513156/SRR15513156.uniq.sorted.bam \
    $BAM_DIR/SRR15513157/SRR15513157.uniq.sorted.bam \
    $BAM_DIR/SRR15513158/SRR15513158.uniq.sorted.bam

##########################
## change to merged bam ##
##########################

samples=(SRR155131{48..58} SRR15513148_49_50_51_52 SRR15513153_54_55_56_57_58)

# 6. Ribo-TISH workflow
echo "### Ribo-TISH workflow ###"
## for regular ribo-seq data
for sample in ${samples[@]};
do
    read -u6
    {
        echo "-- for $sample --"
        [ -d $RiboTISH_DIR/${sample} ] || mkdir -p $RiboTISH_DIR/${sample}
        cd $RiboTISH_DIR/${sample}
        # activate environment
        # conda activate ribotish
        # 6.1 Quality control of riboseq bam data
        echo "-- Quality control of riboseq bam data --"
        # index bam file
        [ -f $BAM_DIR/${sample}/${sample}.uniq.sorted.bam.bai ] || samtools index -@ 20 $BAM_DIR/${sample}/${sample}.uniq.sorted.bam
        [ -f ribotish_qual.pdf ] || ribotish quality \
            -b $BAM_DIR/${sample}/${sample}.uniq.sorted.bam \
            -g $GTF -p 20 \
            -o ./ribotish_qual.txt -f ./ribotish_qual.pdf -r ./ribotish.para.py \
            1>ribotish.quality.log 2>&1

        
        # 6.2 Quality control of riboseq bam data for high-confident de novo genes
        echo "-- Quality control of riboseq bam data for high-confident de novo genes --"
        [ -f ribotish_qual.hc_denovo.pdf ] || ribotish quality \
            -b $BAM_DIR/${sample}/${sample}.uniq.sorted.bam \
            -g $DENOVO_HC_GTF -p 20 \
            -o ./ribotish_qual.hc_denovo.txt -f ./ribotish_qual.hc_denovo.pdf -r ./ribotish.para.hc_denovo.py \
            1>ribotish.quality.hc_denovo.log 2>&1

        # 6.3 Quality control of riboseq bam data for candidate de novo genes
        echo "-- Quality control of riboseq bam data for $candidate --"
        for candidate in ${candidates[@]};
        do
            [ -f ribotish_qual.$candidate.pdf ] || ribotish quality \
                -b $BAM_DIR/${sample}/${sample}.uniq.sorted.bam \
                -g $DENOVO_HC_DIR/gtf/$candidate.gtf -p 20 \
                -o ./ribotish_qual.$candidate.txt -f ./ribotish_qual.$candidate.pdf -r ./ribotish.para.$candidate.py \
                1>ribotish.quality.$candidate.log 2>&1
        done

        # 6.4 Predict ORF/TIS with riboseq bam files (Only test input candidate ORFs)
        echo "-- Predict ORF/TIS with riboseq bam files --"
        [ -f denovo_ORFs.pred.txt ] || ribotish predict \
            -b $BAM_DIR/${sample}/${sample}.uniq.sorted.bam \
            -g $GTF -f $GENOME.fa -i $CandidateORF_DIR/denovo_ORFs.txt \
            --ribopara ./ribotish.para.denovo.py \
            -o denovo_ORFs.pred.txt --tpth 1 --fpth 1 --fspth 1 --fsqth 1 -v -p 20 \
            1>ribotish.predict.log 2>&1
        
        # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
        echo >&6
    }&
done
wait
