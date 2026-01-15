
# de novo ORF以及noncano ORF参考http://10.10.30.30:6767/notebooks/data3/project/sORFs/09-CustomDb/formal_20250821/scripts/generate_evidence_matrix.ipynb
# cano ORF以及lncORF参考http://10.10.30.30:6767/notebooks/data3/project/sORFs/09-CustomDb/formal_20250821/scripts/choose_appropriate_cutoff.ipynb
# 合并@chunfu xiao得到的intergenic结果
# less /home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.bed
# less /home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.fa
# less /home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.aa.fa
# 去掉.fa文件中seq名字中::以及之后的部分，随机抽取小于等于150aa的2500条，并且得到对应的.fa以及.bed文件
# 得到的CDS以及pep都是不包含终止密码子的，同时CDS的bed文件是0-based
source activate biotools
# python seve_list_compare/clean.intergenic.orf.py

##### 一、得到combine的input for ortholog extraction ##### 
# 1、得到CDS bed
# 共有intergenic、lnc、r1m- noncano、noncano、cano物种
# noncano的已经不用拿去做同源基因提取，但是需要提取cano、lnc以及integenic
# 需要得到noncano的CDS bed用来提取vcf
input_dir=../processed/seve_list_compare/
out_dir=$input_dir/input_for_ortholog_extraction
mkdir -p $out_dir
total_orf_gtf=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.gtf
total_orf_gtf_cds_bed=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.CDS.bed
total_orf_gtf_1=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.gtf
total_orf_gtf_cds_bed_1=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.CDS.bed
get_cds_bed(){
	less  $1 | awk '$3=="CDS"'| awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > $2
}
# get_cds_bed $total_orf_gtf $total_orf_gtf_cds_bed
# get_cds_bed $total_orf_gtf_1 $total_orf_gtf_cds_bed_1

suffix=sample.2.5k.lt.150aa

for type in canonical_orfs lnc_orfs;do
    prefix=$input_dir/$type.$suffix
    cut -f 1 $prefix.txt | tail -n +2 -> $prefix.id.txt
    grep -F -f $prefix.id.txt $total_orf_gtf_cds_bed > $prefix.bed
done
for type in non_cano_orfs r1_m_minus_non_cano_orfs;do
    prefix=$input_dir/$type.$suffix
    cut -f 1 $prefix.txt | tail -n +2 -> $prefix.id.txt
    grep -F -f $prefix.id.txt $total_orf_gtf_cds_bed_1 > $prefix.bed
done
# denovo 基因
denovo_prefix=/home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list
grep -F -f $denovo_prefix.id.txt $total_orf_gtf_cds_bed_1 > $denovo_prefix.bed

seqkit seq -n $input_dir/intergenic_orfs.$suffix.aa.fa > $input_dir/intergenic_orfs.$suffix.id.txt
cat $input_dir/canonical_orfs.$suffix.id.txt \
    $input_dir/lnc_orfs.$suffix.id.txt \
    $input_dir/intergenic_orfs.$suffix.id.txt > $input_dir/canonical_orfs.lnc_orfs.intergenic_orfs.id.txt

# 合并
## for kaks同源区域提取
cat $input_dir/intergenic_orfs.$suffix.bed \
    $input_dir/canonical_orfs.$suffix.bed \
    $input_dir/lnc_orfs.$suffix.bed > $out_dir/canonical.intergenic_orfs.lnc_orfs.bed
# 所有（for VCF提取）
cat $input_dir/intergenic_orfs.$suffix.bed \
    $input_dir/lnc_orfs.$suffix.bed \
    $input_dir/canonical_orfs.$suffix.bed \
    $input_dir/non_cano_orfs.$suffix.bed \
    $input_dir/r1_m_minus_non_cano_orfs.$suffix.bed  > $out_dir/intergenic.lnc.noncano.cano.orfs.bed
# 2、得到protein.aa
# 生成FASTA文件
input_dir=../processed/seve_list_compare/
for orfs in lnc_orfs r1_m_minus_non_cano_orfs non_cano_orfs canonical_orfs; do
    cut -f 1,10 $input_dir/$orfs.$suffix.txt | tail -n +2 | seqkit tab2fx > $input_dir/$orfs.$suffix.fa
done
cut -f 1,10 $denovo_prefix.txt | tail -n +2 | seqkit tab2fx > $denovo_prefix.fa
## for kaks同源区域提取
cat $input_dir/intergenic_orfs.$suffix.aa.fa \
    $input_dir/canonical_orfs.$suffix.fa \
    $input_dir/lnc_orfs.$suffix.fa > $input_dir/input_for_pnps/canonical.intergenic_orfs.lnc_orfs.fa

###### 二、整理所有集合的ID、类型 ######
input_dir=../processed/seve_list_compare/
# 合并所有ID并添加类型标签
{
    awk '{print $0"\tIntergenic"}' $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.id.txt
    awk '{print $0"\tlncRNA"}' $input_dir/lnc_orfs.id.txt
    awk '{print $0"\tNon-canonical-ORFs"}' $input_dir/r1_m_minus_non_cano_orfs.sample.2.5k.lt.150aa.id.txt
    awk '{print $0"\tNon-canonical-Protein"}' $input_dir/non_cano_orfs.sample.2.5k.lt.150aa.id.txt
    awk '{print $0"\tCanonical"}' $input_dir/canonical_orfs.sample.2.5k.lt.150aa.id.txt
} | awk 'BEGIN{print "original_id\ttype\tsafe_id"} 
         {
             orig=$1
             type=$2
             safe=orig
             gsub(/[^A-Za-z0-9._-]/, "__", safe)
             gsub(/__+/, "__", safe)
             print orig"\t"type"\t"safe
         }' > $input_dir/all_orfs_with_type.tmp.txt

{
    awk '{print $0"\tDenovo"}' $denovo_prefix.id.txt
} | awk 'BEGIN{print "original_id\ttype\tsafe_id"} 
         {
             orig=$1
             type=$2
             safe=orig
             gsub(/[^A-Za-z0-9._-]/, "_", safe)
             print orig"\t"type"\t"safe
         }' > $input_dir/denovo_orfs_with_type.txt

cat $input_dir/all_orfs_with_type.tmp.txt  $input_dir/denovo_orfs_with_type.txt > $input_dir/all_orfs_with_type.txt
##### 三、得到combine的gtf for snpeff注释 ##### 
out_dir=../processed/seve_list_compare/input_for_pnps
mkdir -p $out_dir
suffix=sample.2.5k.lt.150aa
total_orf_gtf=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.gtf
total_orf_gtf_1=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.gtf
for type in canonical_orfs lnc_orfs;do
    prefix=$input_dir/$type.$suffix
    cut -f 1 $prefix.txt | tail -n +2 -> $prefix.id.txt
    grep -F -f $prefix.id.txt $total_orf_gtf > $prefix.gtf
done
for type in non_cano_orfs r1_m_minus_non_cano_orfs;do
    prefix=$input_dir/$type.$suffix
    cut -f 1 $prefix.txt | tail -n +2 -> $prefix.id.txt
    grep -F -f $prefix.id.txt $total_orf_gtf_1 > $prefix.gtf
done
grep -F -f $denovo_prefix.id.txt $total_orf_gtf_1 > $denovo_prefix.gtf
## 添加intergenic的gtf
awk 'BEGIN{OFS="\t"} {
    if ($6 == "+") {
        start = $2 + 1
        end = $3 + 3  # 包含终止密码子
    } else {
        start = $2 - 2  # 包含终止密码子
        end = $3
    }
    
    attr = "gene_id \"" $4 "\"; transcript_id \"" $4 "\";"
    
    # transcript
    print $1, "intergenic_ORF", "transcript", start, end, ".", $6, ".", attr
    # exon
    print $1, "intergenic_ORF", "exon", start, end, ".", $6, ".", attr
    # CDS
    print $1, "intergenic_ORF", "CDS", start, end, ".", $6, "0", attr
}' $input_dir/intergenic_orfs.$suffix.bed > $input_dir/intergenic_orfs.$suffix.gtf

>$out_dir/all_orfs.gtf
for type in canonical_orfs lnc_orfs non_cano_orfs r1_m_minus_non_cano_orfs intergenic_orfs;do
    prefix=$input_dir/$type.$suffix
    ls $prefix.gtf
    cat $prefix.gtf >> $out_dir/all_orfs.gtf
done
# less $out_dir/all_orfs.gtf |awk '$3=="transcript"'|wc -l

##### 四、得到combine的cds序列以及protein序列 ##### 
out_dir=../processed/seve_list_compare/input_for_pnps
# cp $input_dir/intergenic_orfs.$suffix.fa $input_dir/intergenic_orfs.$suffix.cds.fa 
ln -s $(realpath $input_dir/intergenic_orfs.$suffix.aa.fa) $input_dir/intergenic_orfs.$suffix.fa
# protein fa
>$out_dir/all_orfs.fa
for type in canonical_orfs lnc_orfs non_cano_orfs r1_m_minus_non_cano_orfs intergenic_orfs;do
    prefix=$input_dir/$type.$suffix
    # ls $prefix.fa
    # head -n2 $prefix.fa
    cat $prefix.fa >> $out_dir/all_orfs.fa
done

# CDS fa
>$out_dir/all_orfs.cds.tmp.fa
for type in lnc_orfs r1_m_minus_non_cano_orfs non_cano_orfs canonical_orfs; do
    prefix=$input_dir/$type.$suffix
    cut -f 1,28 $prefix.txt | tail -n +2 | seqkit tab2fx > $prefix.cds.fa
    ls $prefix.cds.fa
    cat $prefix.cds.fa >> $out_dir/all_orfs.cds.tmp.fa
done
cut -f 1,28 $denovo_prefix.txt | tail -n +2 | seqkit tab2fx > $denovo_prefix.cds.fa
cat $out_dir/all_orfs.cds.tmp.fa $input_dir/intergenic_orfs.$suffix.fa $denovo_prefix.cds.fa > $out_dir/all_orfs.cds.fa

# rm -rf $out_dir/all_orfs.cds.tmp.fa