###### S1 ######
# 根据sorf id补充信息
Rscript 1.Uni.annotate.basic.v1.20250528.R "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/cano_sep_orf_id.txt" \
"FALSE" "FALSE" "$PWD/output/S1/" "cano_sep_tab_info.txt"

# 首先去掉在galaxy上导出出现一些小bug的行，例如两个文件合并的时候都有列名，以及有一些Protein列为空（因为有些肽段没有比对到其他蛋白质上，但是is.unique列为false）
awk -F'\t' -v OFS='\t' '$35 != "\"\""' /home/user/data3/lit/project/sORFs/02-Mass-spec/transferred_from_galaxy/20250528/psm_sep_all.txt |grep -v Protein > $PWD/output/S1/psm_sep_all.txt
cat $PWD/output/S1/psm_sep_all.txt|awk -v FS='\t' -v OFS='\t' '{print $35}'|grep -v sp|sort -u > $PWD/output/S1/noncano.sep.ms.specific.nonspecific.tmp.lst
sed '1i ORF_id_trans' $PWD/output/S1/noncano.sep.ms.specific.nonspecific.tmp.lst > $PWD/output/S1/noncano.sep.ms.specific.nonspecific.lst
rm -rf $PWD/output/S1/noncano.sep.ms.specific.nonspecific.tmp.lst
# 29268
wc -l $PWD/output/S1/noncano.sep.ms.specific.nonspecific.lst
nohup Rscript 1.Uni.annotate.basic.v1.20250528.R "$PWD/output/S1/noncano.sep.ms.specific.nonspecific.lst" \
"TRUE" "TRUE" "$PWD/output/S1/" "uncano_sep_tab_info.txt" &

###### S2 ######
# 1. 根据bam推测文库类型
bash 2.Uni.libType.v1.20250529.sh /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227/Aligned.sortedByCoord.out.bam.lst output/S2/Ribo-seq
# 2. 整理bam lst
2.1*R
# 3-5
for i in rna ribo;do
# for i in ribo;do
# 3. featureCounts
output_path=output/S2/fc_output_${i}_seq
mkdir -p $output_path && cd $output_path
gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
fc_script=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/2.2.Uni.fc.20250530.sh
lst=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/output/S2/total_${i}_bam_lst.txt
bash $fc_script -i $lst -g $gtf -o ${i}_seq &> fc.log
# cd -
# 4. libsize
libsize_script=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250528/Uni.calcu.libsize.v1.20250530.sh
cut -f 1 $lst | tail -n +2 > output/S2/total_${i}_bam_lst.1.txt
bash $libsize_script output/S2/total_${i}_bam_lst.1.txt $output_path
# 5. rpkm
Rscript 2.3.Uni.get.rpkm.20250530.R $PWD/$output_path/${i}_seq_combined.txt \
 $PWD/$output_path/libsize.txt $lst
done

# 2.Uni.libType.v1.20250529.sh可以增强SingleEnd的输出，目前是手动赋值
# 6. 整理样本信息，查看样本之间的相关性
2.4*
2.5*
# 7. 合并相同研究相同发育时期的样本，探究RNA-seq以及Ribo-seq之间的相关性
2.6*