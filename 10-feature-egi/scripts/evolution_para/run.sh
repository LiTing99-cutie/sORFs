# 0.1 导出自己注释的ORF，以便和经典的注释合并
out_dir=../processed/seve_list_compare/
out_dir_1=../processed/seve_list_compare/input_for_pnps/
all_r1_m1_orf_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf
less /home/user/data3/lit/project/sORFs/11-denovo-list/processed/non_cano_orfs_r1m1.txt |cut -f 1|tail -n  +2|sort -u > $out_dir/non_cano_orfs_r1m1.id.txt
grep -F -f <(awk '{print "gene_id \""$1"\""}' $out_dir/non_cano_orfs_r1m1.id.txt) \
    $all_r1_m1_orf_gtf \
    > $out_dir/non_cano_orfs_r1m1.orfs.gtf
# 进一步和lnc以及intergenic ORF合并
cat $out_dir/non_cano_orfs_r1m1.orfs.gtf $out_dir/intergenic_orfs.sample.2.5k.lt.150aa.gtf \
 $out_dir/lnc_orfs.gtf > $out_dir_1/custom.gtf

# 计算de novo ORF以及canonical ORF以及lncRNA ORF的ka/ks以及pn/ps
## kaks需要提取同源的蛋白质以及CDS序列【de novo基因已经有】；需要检查人中的蛋白质序列
## 是否和*ortholog.pep.fa中的匹配；以及需要把首个氨基酸替换成M再进行匹配的比较
## pnps可以直接根据基因组的坐标计算

##### pnps #####
# 1.1 注释vcf
nohup bash pnps/1.annotate_vcf.sh > ../log/annotate_vcf.log 2>&1 &
# nohup bash evolution_para/1.annotate_vcf.test.sh > ../log/annotate_vcf.test.log 2>&1 &
# tail -f ../log/annotate_vcf.log
# 1.2 提取目标区域的vcf
# vcf=/home/user/data3/lit/project/sORFs/11-denovo-list/public_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_AllVar_QC_metrics.vcf.gz
# chr1    99335318        99335387        99309590-99464377_203   .       +
nohup bash pnps/2.extract.vcf.sh &> ../log/extract.vcf.log &
nohup bash pnps/3.annotate.vcf.sh &> ../log/annotate.vcf.log &

##### kaks #####
# # 1.1 单个case测试
# ## 需要修改paraAT的脚本
# denovo_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251029_denovo_check/processed
# orf=PB.29415.11__chr22__-__201__3603__2541__2736__dORF__CTG
# bash ./kaks/kaks.bin.sh \
#     -p $denovo_dir/check_denovo/chr22/orfs/$orf.ortholog.pep.fa \
#     -c $denovo_dir/check_denovo/chr22/orfs/$orf.ortholog.nucl.fa \
#     -o ../processed/evolution_para/kaks/$orf \
#     -m human \
#     -t 10
# # 1.2 批量运行和合并
# ## 1.2.1 小批量试运行
# find $denovo_dir -name "*ortholog.pep.fa" | head -n50 | while read pep; do
#     nucl="${pep%.pep.fa}.nucl.fa"
#     cp "$pep" tmp/
#     cp "$nucl" tmp/
# done
# nohup bash ./kaks/kaks.parallel.sh $PWD/tmp ../processed/evolution_para/kaks &> ../log/kaks.parallel.log &

## 1.2.2 正式运行
nohup bash kaks/cp.noncano.ortholog.res.sh &> ../log/cp.noncano.ortholog.res.log &
nohup bash kaks/cp.ortholog.res.translate.sh &> ../log/cp.ortholog.res.translate.log &
nohup bash kaks/kaks.run.2025114.sh &> ../log/kaks.run.20251114.log &
# 监控：
# ls /home/user/data3/lit/project/sORFs/11-denovo-list/processed/evolution_para/kaks|wc -l