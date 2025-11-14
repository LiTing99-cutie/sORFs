
# de novo ORF以及noncano ORF参考http://10.10.30.30:6767/notebooks/data3/project/sORFs/09-CustomDb/formal_20250821/scripts/generate_evidence_matrix.ipynb
# cano ORF以及lncORF参考http://10.10.30.30:6767/notebooks/data3/project/sORFs/09-CustomDb/formal_20250821/scripts/choose_appropriate_cutoff.ipynb
# 合并@chunfu xiao得到的intergenic结果
# less /home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.bed
# less /home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.fa
# less /home/user/data3/rbase/denovo_tumor/denovo_genes/all_ORFs/intergenic_ORFs/intergenic_ORFs.gt30.shuffle.aa.fa
# 去掉.fa文件中seq名字中::以及之后的部分，随机抽取小于等于150aa的2500条，并且得到对应的.fa以及.bed文件
# 得到的CDS以及pep都是不包含终止密码子的，同时CDS的bed文件是0-based
source activate biotools
python seve_list_compare/clean.intergenic.orf.py

##### 一、得到combine的input for ortholog extraction ##### 
# 1、得到CDS bed
# noncano的已经不用拿去做同源基因提取，但是需要提取sample的经典ORF以及lncORF
# 需要得到noncano的CDS bed用来提取vcf
input_dir=../processed/seve_list_compare/
out_dir=../processed/seve_list_compare/input_for_ortholog_extraction
mkdir -p $out_dir
total_orf_gtf=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.gtf
total_orf_gtf_cds_bed=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v3_20251103/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.CDS.bed
total_orf_gtf_1=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.gtf
total_orf_gtf_cds_bed_1=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/candidateORF.filtered.rmdup.renamed.pep.extended.1.txt/sub.CDS.bed
get_cds_bed(){
	less  $1 | awk '$3=="CDS"'| awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > $2
}
get_cds_bed $total_orf_gtf $total_orf_gtf_cds_bed
get_cds_bed $total_orf_gtf_1 $total_orf_gtf_cds_bed_1

for prefix in $input_dir/canonical_orfs.sample.2.5k.lt.150aa $input_dir/lnc_orfs;do
    cut -f 1 $prefix.txt | tail -n +2 -> $prefix.id.txt
    grep -F -f $prefix.id.txt $total_orf_gtf_cds_bed > $prefix.bed
done
for prefix in $input_dir/non_cano_orfs.sample.2.5k.lt.150aa;do
    cut -f 1 $prefix.txt | tail -n +2 -> $prefix.id.txt
    grep -F -f $prefix.id.txt $total_orf_gtf_cds_bed_1 > $prefix.bed
done
seqkit seq -n $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.aa.fa > $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.id.txt
cat $input_dir/canonical_orfs.sample.2.5k.lt.150aa.id.txt \
    $input_dir/lnc_orfs.id.txt \
    $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.id.txt > $input_dir/canonical_orfs.lnc_orfs.intergenic_orfs.id.txt

# 合并
cat $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.bed \
    $input_dir/canonical_orfs.sample.2.5k.lt.150aa.bed \
    $input_dir/lnc_orfs.bed > $out_dir/canonical.intergenic_orfs.lnc_orfs.bed
cat $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.bed \
    $input_dir/canonical_orfs.sample.2.5k.lt.150aa.bed \
    $input_dir/lnc_orfs.bed \
    $input_dir/non_cano_orfs.sample.2.5k.lt.150aa.bed > $out_dir/noncano.canonical.intergenic_orfs.lnc_orfs.bed
# 2、得到protein.aa
# 生成FASTA文件
for orfs in canonical_orfs; do
    cut -f 1,10 $input_dir/$orfs.sample.2.5k.lt.150aa.txt | tail -n +2 | seqkit tab2fx > $input_dir/$orfs.sample.2.5k.lt.150aa.fa
done
cat $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.aa.fa \
    $input_dir/canonical_orfs.sample.2.5k.lt.150aa.fa \
    $input_dir/lnc_orfs.fa > $out_dir/canonical.intergenic_orfs.lnc_orfs.fa

##### 一、得到combine的gtf for snpeff注释 ##### 
out_dir=../processed/seve_list_compare/input_for_pnps
mkdir $out_dir
# 经典的和非经典的都不需要，只需要添加intergenic以及lncRNA的gtf
## 添加lnc_orfs的gtf
prefix=../processed/seve_list_compare/lnc_orfs
# 提取ID列表
cut -f1 $prefix.txt | tail -n +2 > ${prefix}_ids.tmp
# 一次性grep所有ID
grep -F -f ${prefix}_ids.tmp $total_orf_gtf > $prefix.gtf
# 清理临时文件
rm -rf ${prefix}_ids.tmp
# less $prefix.gtf |awk '$3=="transcript"'|wc -l

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
}' $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.bed > $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.gtf
cat $prefix.gtf $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.gtf > $out_dir/intergenic_orfs.lnc_orfs.gtf
# less $input_dir/intergenic_orfs.sample.2.5k.lt.150aa.gtf |awk '$3=="transcript"'|wc -l
