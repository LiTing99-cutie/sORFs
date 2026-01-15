################################################
#File Name: run.pnps.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Fri 01 Mar 2024 11:54:58 AM CST
################################################

#!/bin/sh
##并发运行脚本，并控制并发数
# 设置并发的进程数
thread_num=10
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
VCF_FILE=/home/user/data/rbase/Genome_data/human/GTExV8_WGS/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites_VEP_annot.vcf.gz
TOOL_DIR=/home/user/data3/rbase/opt/snpEff
LIST_DIR=/home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs
GTF=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/Homo_sapiens.GRCh38.gencode.v43.annotation.gtf
DENOVO_GPE=/home/user/data3/rbase/denovo_tumor/denovo_genes/100_denovo_gene_list/100_denovo_genes.compiled.gpe
DENOVO_LIST=/home/user/data3/rbase/denovo_tumor/denovo_genes/100_denovo_gene_list/100_denovo_genes.compiled.txt
ALL_GTF=$TOOL_DIR/data/hg38/genes.gtf
JAVA_DIR=/usr/lib/jvm/java-11-openjdk-11.0.22.0.7-1.el7_9.x86_64/bin

###################################################################################
## Reconstruct gtf file containing noncoding and de novo genes (ORF information) ##
###################################################################################
echo "## Reconstruct gtf file containing noncoding and de novo genes ##"

# genepred file (gt30)
cut -f 1-10 $LIST_DIR/pc_nc_ORFs.gt30.genepred.txt > $LIST_DIR/all_ORFs.gt30.genepred.txt
# add intron, intergenic regions and de novo genes
awk 'BEGIN{OFS=FS="\t"}{print $4, $1, $6, $2, $3, $2, $3, 1, $2",", $3","}' $LIST_DIR/intron_ORFs/intron_ORFs.gt30.shuffle.bed \
    >> $LIST_DIR/all_ORFs.gt30.genepred.txt
awk 'BEGIN{OFS=FS="\t"}{print $4, $1, $6, $2, $3, $2, $3, 1, $2",", $3","}' $LIST_DIR/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.bed \
    >> $LIST_DIR/all_ORFs.gt30.genepred.txt
awk 'BEGIN{OFS=FS="\t"}{print $4, $1, $6, $2, $3, $2, $3, 1, $2",", $3","}' $LIST_DIR/intergenic_random_seq/intergenic_random.bed \
    >> $LIST_DIR/all_ORFs.gt30.genepred.txt
cat $DENOVO_GPE >> $LIST_DIR/all_ORFs.gt30.genepred.txt

# generate gtf file
cat $LIST_DIR/all_ORFs.gt30.genepred.txt | genePredToGtf file stdin $LIST_DIR/all_ORFs.gt30.gtf
cp $LIST_DIR/all_ORFs.gt30.gtf $ALL_GTF

######################################
## Run snpEff to annotate vcf files ##
######################################
# generate bed file
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4}' $LIST_DIR/all_ORFs.gt30.ORF_no_intron.bed > $WORK_DIR/all_ORFs.gt30.ORF_no_intron.bed4
cd $TOOL_DIR
## filter focused variants
echo "## Filter off focused variants ##"
[ -f $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites_VEP_annot.all_ORFs.gt30.vcf.gz ] || \
bcftools filter -Oz --threads 20 $VCF_FILE -R $WORK_DIR/all_ORFs.gt30.ORF_no_intron.bed4 \
    -o $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites_VEP_annot.all_ORFs.gt30.vcf.gz
## index
[ -f $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites_VEP_annot.all_ORFs.gt30.vcf.gz.csi ] || \
    bcftools index --threads 10 $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites_VEP_annot.all_ORFs.gt30.vcf.gz

# run snpEff to build database
echo "## Run snpEff to build database ##"
## -noCheckCds: Skip CDS sequences check.
## -noCheckProtein: Skip Protein sequences check.
## -noCheckCds -noCheckProtein: When using BOTH command line options, SnpEff to save the new database without checking it (i.e. neither CDS nor protein sequences are checked).
[ -f $TOOL_DIR/data/hg38/snpEffectPredictor.bin ] || $JAVA_DIR/java -Xmx4g -jar $TOOL_DIR/snpEff.jar build \
    -c $TOOL_DIR/snpEff.config -gtf22 -v -noCheckCds -noCheckProtein hg38

# run snpEff to annotate vcf file
echo "## Run snpEff to annotate vcf file ##"
[ -f $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz ] || \
    $JAVA_DIR/java -Xmx10G -jar $TOOL_DIR/snpEff.jar eff -c $TOOL_DIR/snpEff.config -v \
        hg38 $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites_VEP_annot.all_ORFs.gt30.vcf.gz \
        > $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf -csvStats GTEx.csv -stats GTEx.html

# bgzip
[ -f $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz ] || \
    bcftools view $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf \
        -Oz -o $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz
[ -f $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf ] && \
    rm $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf
# index 
[ -f $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz.csi ] || \
    bcftools index --threads 20 $TOOL_DIR/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_sites.snpEff.all_ORFs.gt30.vcf.gz

##########################################################
## calculate expected missense sites and synonymy sites ##
##########################################################

# python cal_n_s_sites.py --output_file ./all.gt30.orf_seq.n_s_sites &

################################
## extract syn and nonsyn VCF ##
################################
mkdir -p $WORK_DIR/pnps/synonymous_mutation
mkdir -p $WORK_DIR/pnps/missense_variant
mkdir -p $WORK_DIR/bed
mkdir -p $WORK_DIR/vcf_files

## for each ORF
echo -e "ORF ID\tNum_Syn_Var\tFreq_Syn_Var\tNum_Mis_Var\tFreq_Mis_Var" > $WORK_DIR/pnps/all_ORFs.gt30.syn_mis.num.txt
echo -e "ORF ID\tNum_Syn_Var\tFreq_Syn_Var\tNum_Mis_Var\tFreq_Mis_Var" > $WORK_DIR/pnps/all_ORFs.gt30.syn_mis.num.pos.txt
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
        
        # synonymous mutation
        # DB: dbSNP Membership
        # POSITIVE_TRAIN_SITE (VQSR): This variant was used to build the positive training set of good variants.
        # NEGATIVE_TRAIN_SITE (VQSR): This variant was used to build the negative training set of bad variants
        # | grep ";DB;" | grep "POSITIVE_TRAIN_SITE"
        zcat $WORK_DIR/vcf_files/tmp.$orf.vcf.gz | grep -v "^#" | grep ";DB;" | grep "synonymous_variant" | \
            awk 'BEGIN{OFS=FS="\t"}
            {
                vep_i=match($8,"CSQ=")
                snpeff_i=match($8,"ANN=");
                split(substr($8,snpeff_i),snpeff,"|,")
                for(i=1;i<=length(snpeff);i++)
                {
                    print $1,$2,substr($8,0,vep_i),snpeff[i]
                }
            }' | grep $orf | grep "synonymous_variant" \
            > $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.tmp
        ## POSITIVE_TRAIN_SITE
        grep "POSITIVE_TRAIN_SITE" $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.tmp\
            > $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.pos.tmp
        # count
        num_syn_mu=`cat $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.tmp | wc -l`
        freq_syn_mu=`grep -oE ";AF=0.[0-9]*;" $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.tmp | sed "s@AF=@@g" | sed "s@;@@g" | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
        num_syn_mu_pos=`cat $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.pos.tmp | wc -l`
        freq_syn_mu_pos=`grep -oE ";AF=0.[0-9]*;" $WORK_DIR/pnps/synonymous_mutation/synonymous_mutation-$orf.pos.tmp | sed "s@AF=@@g" | sed "s@;@@g" | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
        
        # missense_variant mutation
        zcat $WORK_DIR/vcf_files/tmp.$orf.vcf.gz | grep -v "^#" | grep ";DB;" | grep "missense_variant" | \
            awk 'BEGIN{OFS=FS="\t"}
            {
                vep_i=match($8,"CSQ=")
                snpeff_i=match($8,"ANN=");
                split(substr($8,snpeff_i),snpeff,"|,")
                for(i=1;i<=length(snpeff);i++)
                {
                    print $1,$2,substr($8,0,vep_i),snpeff[i]
                }
            }' | grep $orf | grep "missense_variant" \
            > $WORK_DIR/pnps/missense_variant/missense_variant-$orf.tmp
        ## POSITIVE_TRAIN_SITE
        grep "POSITIVE_TRAIN_SITE" $WORK_DIR/pnps/missense_variant/missense_variant-$orf.tmp \
            > $WORK_DIR/pnps/missense_variant/missense_variant-$orf.pos.tmp
        # count
        num_mis_mu=`cat $WORK_DIR/pnps/missense_variant/missense_variant-$orf.tmp | wc -l`
        freq_mis_mu=`grep -oE ";AF=0.[0-9]*;" $WORK_DIR/pnps/missense_variant/missense_variant-$orf.tmp | sed "s@AF=@@g" | sed "s@;@@g" | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`
        num_mis_mu_pos=`cat $WORK_DIR/pnps/missense_variant/missense_variant-$orf.pos.tmp | wc -l`
        freq_mis_mu_pos=`grep -oE ";AF=0.[0-9]*;" $WORK_DIR/pnps/missense_variant/missense_variant-$orf.pos.tmp | sed "s@AF=@@g" | sed "s@;@@g" | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'`

        # summary
        echo -e "$orf\t$num_syn_mu\t$freq_syn_mu\t$num_mis_mu\t$freq_mis_mu" >> $WORK_DIR/pnps/all_ORFs.gt30.syn_mis.num.txt
        echo -e "$orf\t$num_syn_mu_pos\t$freq_syn_mu_pos\t$num_mis_mu_pos\t$freq_mis_mu_pos" >> $WORK_DIR/pnps/all_ORFs.gt30.syn_mis.num.pos.txt

        echo >&6
    }&
done
wait