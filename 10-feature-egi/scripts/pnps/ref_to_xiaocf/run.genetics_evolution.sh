################################################
#File Name: run.avg_summary.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Wed 24 Jan 2024 01:42:48 PM CST
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

WORK_DIR=/home/user/data3/rbase/denovo_tumor/denovo_genes/genetics
TOOL_DIR=/home/user/data3/rbase/opt/snpEff
VCF_DIR=/home/user/data/rbase/Genome_data/human/1000_genomes_project
LIST_DIR=/home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs

awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4}' $LIST_DIR/all_ORFs.gt30.ORF_no_intron.bed > $WORK_DIR/all_ORFs.gt30.ORF_no_intron.bed4
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4}' $LIST_DIR/all_ORFs.gt100.ORF_no_intron.bed > $WORK_DIR/all_ORFs.gt100.ORF_no_intron.bed4
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4}' $LIST_DIR/all_ORFs.gt300.ORF_no_intron.bed > $WORK_DIR/all_ORFs.gt300.ORF_no_intron.bed4
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4}' $LIST_DIR/all_ORFs.gt600.ORF_no_intron.bed > $WORK_DIR/all_ORFs.gt600.ORF_no_intron.bed4


    #########################
    ## phyloP for all ORFs ##
    #########################

    # cactus241way
    [ -f $WORK_DIR/cactus241way/all_ORFs.gt30.ORF_no_intron.tab ] || bigWigAverageOverBed $WORK_DIR/cactus241way/hg38.cactus241way.phyloP.bw \
        $WORK_DIR/all_ORFs.gt30.ORF_no_intron.bed4 $WORK_DIR/cactus241way/all_ORFs.gt30.ORF_no_intron.tab

    # phyloP100way
    [ -f $WORK_DIR/phyloP100way/all_ORFs.gt30.ORF_no_intron.tab ] || bigWigAverageOverBed $WORK_DIR/phyloP100way/hg38.phyloP100way.bw \
        $WORK_DIR/all_ORFs.gt30.ORF_no_intron.bed4 $WORK_DIR/phyloP100way/all_ORFs.gt30.ORF_no_intron.tab

if [[ 1 > 2 ]];then
    #######################################
    ## Nucleotide diversity for all ORFs ##
    #######################################

    # pi for sites
    [ -f $WORK_DIR/pi/nucleotide_diversity.sites.pi ] || \
        zcat $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz | \
        vcftools --vcf - --site-pi --bed $WORK_DIR/all_ORFs.gt30.ORF_no_intron.bed4 --out $WORK_DIR/pi/nucleotide_diversity

    # pi for ORFs
    mkdir -p $WORK_DIR/vcf_files
    mkdir -p $WORK_DIR/bed
    mkdir -p $WORK_DIR/pi/tmp
    mkdir -p $WORK_DIR/pi/logs
    ## for each ORF
    cut -f1 $LIST_DIR/all_ORFs.gt30.ORF_no_intron.saf | uniq | while read orf
    do
        read -u6
        {
            # extract bed of the orf
            grep $orf $LIST_DIR/all_ORFs.gt30.ORF_no_intron.saf | awk 'BEGIN{OFS=FS="\t"}{print $2,$3,$4,$1}' > $WORK_DIR/bed/tmp.$orf.bed4
            echo " -- for $orf -- "
        
            # extract vcf annotations of the orf
            [ -f $WORK_DIR/vcf_files/tmp.$orf.vcf.gz ] || \
            bcftools filter -Oz --threads 20 $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz \
                -R $WORK_DIR/bed/tmp.$orf.bed4 > $WORK_DIR/vcf_files/tmp.$orf.vcf.gz

            # calculate pi
            [ -f $WORK_DIR/pi/tmp/nucleotide_diversity.$orf.windowed.pi ] || \
            vcftools --vcf $WORK_DIR/vcf_files/tmp.$orf.vcf.gz \
                --window-pi 500000 --out $WORK_DIR/pi/tmp/nucleotide_diversity.$orf \
                1>$WORK_DIR/pi/logs/vcftools.$orf.log 2>&1

            echo >&6
        }&
    done
    wait

    ## merge pi
    echo -e "ORF\tCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" > $WORK_DIR/pi/nucleotide_diversity.orfs.pi
    cut -f1 $LIST_DIR/all_ORFs.gt30.ORF_no_intron.saf | uniq | while read orf
    do
        if [ -f $WORK_DIR/pi/tmp/nucleotide_diversity.$orf.windowed.pi ];then
            awk -v orf=$orf 'BEGIN{OFS=FS="\t"}NR==2{print orf,$0}' $WORK_DIR/pi/tmp/nucleotide_diversity.$orf.windowed.pi \
                >> $WORK_DIR/pi/nucleotide_diversity.orfs.pi
        else
            echo " ---- No nucleotide diversity result for $orf ---- "
        fi
    done
fi

    ########################
    ## pN/pS for all ORFs ##
    ########################

    # generate homologs.txt
    cut -f1 $LIST_DIR/all_ORFs.gt30.ORF_no_intron.saf | uniq | awk 'BEGIN{OFS=FS="\t"}{print $1"_1",$1"_2"}' \
        > $WORK_DIR/pnps/all_ORFs.gt30.homologs.txt
    ## replace special words
    sed -i 's/|/#/g' $WORK_DIR/pnps/all_ORFs.gt30.homologs.txt
    sed -i 's/+/--/g' $WORK_DIR/pnps/all_ORFs.gt30.homologs.txt

    # generate cds.fa
    awk 'BEGIN{OFS=FS="\t"}FNR!=1{
            split($1,id," ")
            print ">"id[1]"_1";
            print $2;
            print ">"id[1]"_2";
            print $2;
        }' \
        $LIST_DIR/coding_ORFs/pc_ORFs.shuffle.fa.txt \
        $LIST_DIR/lncRNA/lncRNA_ORFs.gt30.shuffle.fa.txt \
        $LIST_DIR/uORFs/uORFs.gt30.shuffle.fa.txt \
        $LIST_DIR/dORFs/dORFs.gt30.shuffle.fa.txt \
        $LIST_DIR/intron_ORFs/intron_ORFs.gt30.shuffle.fa.txt \
        $LIST_DIR/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.fa.txt \
        $LIST_DIR/intergenic_random_seq/intergenic_random.fa.txt \
        /home/user/data3/rbase/denovo_tumor/denovo_genes/100_denovo_gene_list/100_denovo_genes.cds.txt \
        > $WORK_DIR/pnps/all_ORFs.gt30.cds.fa
    ## replace special words
    sed -i 's/|/#/g' $WORK_DIR/pnps/all_ORFs.gt30.cds.fa
    sed -i 's/+/--/g' $WORK_DIR/pnps/all_ORFs.gt30.cds.fa

    # generate pep.fa
    faTrans $WORK_DIR/pnps/all_ORFs.gt30.cds.fa $WORK_DIR/pnps/all_ORFs.gt30.pep.fa

    # align by ParaAT.pl
    rm -rf $WORK_DIR/pnps/results
    ParaAT.pl -h $WORK_DIR/pnps/all_ORFs.gt30.homologs.txt -n $WORK_DIR/pnps/all_ORFs.gt30.cds.fa \
        -a $WORK_DIR/pnps/all_ORFs.gt30.pep.fa -p $WORK_DIR/pnps/proc -m mafft -f axt -g -k -o $WORK_DIR/pnps/results

    # calculate S-sites N-sites by KaKs_Calculator
    if [ ! -f $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt ]; then
        echo -e "Sequence\tMethod\tLength\tS_Sites\tN_Sites" > $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt.ng
        rm -rf $WORK_DIR/pnps/ng && mkdir -p $WORK_DIR/pnps/ng
        cat $WORK_DIR/pnps/all_ORFs.gt30.homologs.txt | while read ids
        do
            id1=`echo $ids | cut -f1 -d " "`
            id2=`echo $ids | cut -f2 -d " "`
            echo " -- KaKs_Calculator $id1 $id2 -- "
            KaKs_Calculator -i $WORK_DIR/pnps/results/$id1-$id2.cds_aln.axt -o $WORK_DIR/pnps/ng/$id1-$id2.cds_aln.axt.ng -m GMYN
            tail -n +2 $WORK_DIR/pnps/ng/$id1-$id2.cds_aln.axt.ng | cut -f1,2,7,8,9 >> $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt.ng
        done
    fi

    ## restore special words
    sed -i 's/#/|/g' $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt.ng
    sed -i 's/--/+/g' $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt.ng
    awk 'BEGIN{OFS=FS="\t"}{split($1,A,"_1-"); print A[1],$2,$3,$4,$5}' $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt.ng \
        > $WORK_DIR/pnps/all_ORFs.gt30.cds_aln.axt.txt