# 0.1 导出自己注释的ORF，以便和经典的注释合并
out_dir=../processed/seve_list_compare/
out_dir_1=../processed/seve_list_compare/input_for_pnps/
all_r1_m1_orf_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf
less /home/user/data3/lit/project/sORFs/11-denovo-list/processed/non_cano_orfs_r1m1.txt |cut -f 1|tail -n  +2|sort -u > $out_dir/non_cano_orfs_r1m1.id.txt
grep -F -f <(awk '{print "gene_id \""$1"\""}' $out_dir/non_cano_orfs_r1m1.id.txt) \
    $all_r1_m1_orf_gtf \
    > $out_dir/non_cano_orfs_r1m1.orfs.gtf
# 进一步和lnc以及intergenic ORF合并【gencode中的gtf天然包含所有的经典ORF】
cat $out_dir/non_cano_orfs_r1m1.orfs.gtf $out_dir/intergenic_orfs.sample.2.5k.lt.150aa.gtf \
 $out_dir/lnc_orfs.gtf > $out_dir_1/custom.gtf

# 计算de novo ORF以及canonical ORF以及lncRNA ORF的ka/ks以及pn/ps
## kaks需要提取同源的蛋白质以及CDS序列【de novo基因已经有】；需要检查人中的蛋白质序列
## 是否和*ortholog.pep.fa中的匹配；以及需要把首个氨基酸替换成M再进行匹配的比较
## pnps可以直接根据基因组的坐标计算

##### pnps #####
# 1.1 注释vcf
nohup bash pnps/1.build.annotation.sh > ../log/1.build.annotation.log 2>&1 &
# nohup bash evolution_para/1.annotate_vcf.test.sh > ../log/annotate_vcf.test.log 2>&1 &
# tail -f ../log/annotate_vcf.log
# 1.2 提取目标区域的vcf
# vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
# chr1    99335318        99335387        99309590-99464377_203   .       +
nohup bash pnps/2.extract.vcf.sh &> ../log/extract.vcf.log &
nohup bash pnps/3.annotate.vcf.sh &> ../log/annotate.vcf.log &

# 1.3 计算目标区域的同义位点和非同义位点数量
# vcf=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/pnps/orf_regions/hg38_custom.eff.vcf
# orf_bed=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_pnps/intergenic.lnc.noncano.cano.orfs.bed
bash pnps/4.intersect.orf.bed.vcf.sh