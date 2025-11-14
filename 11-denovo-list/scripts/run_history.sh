# 0.0 导出de novo ORF所在的基因的gtf，查看目的基因在目的细胞系的RNA-seq数据的表达量
# ../processed/denovo_list.txt以及non_cano_orfs_r1m1.txt都是从http://10.10.30.30:6767/notebooks/data3/project/sORFs/09-CustomDb/formal_20250821/scripts/generate_evidence_matrix.ipynb中得到
less ../processed/denovo_list.txt |cut -f 13|tail -n  +2|sort -u > ../processed/denovo_orfs.genes.txt

iso_seq_gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf

grep -f <(sed 's/^/gene_id "/; s/$/"/' ../processed/denovo_orfs.genes.txt) $iso_seq_gtf > ../processed/denovo_orfs.genes.gtf

grep -oP 'gene_id "\K[^"]+' ../processed/denovo_orfs.genes.gtf | sort -u | wc -l

## 需要导出所有的转录本
map=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/representative.tsv
cut -f1 ../processed/denovo_list.txt | tail -n +2 > ../processed/denovo_list.id.txt
# awk 'NR==FNR {ids[$1]=1; next} $2 in ids' \
#     ../processed/denovo_list.id.txt \
#     $map > ../processed/denovo_list.deCollape.orf_id.txt
iso_seq_gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
python extract_transcript_gtf.py \
    ../processed/denovo_list.id.txt \
    $map \
    $iso_seq_gtf \
    ../processed/
# 得到../processed/denovo_orfs.transcripts.complete.gtf
grep -w "gene_id" ../processed/denovo_orfs.transcripts.complete.gtf| awk -F 'gene_id "' '{print $2}' | cut -d '"' -f1 | sort -u | wc -l

# 0.1 导出自己注释的ORF，以便和经典的注释合并
all_r1_m1_orf_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf
less ../processed/non_cano_orfs_r1m1.txt |cut -f 1|tail -n  +2|sort -u > ../processed/non_cano_orfs_r1m1.id.txt
grep -F -f <(awk '{print "gene_id \""$1"\""}' ../processed/non_cano_orfs_r1m1.id.txt) \
    $all_r1_m1_orf_gtf \
    > ../non_cano_orfs_r1m1.orfs.gtf
# 计算de novo ORF以及canonical ORF以及lncRNA ORF的ka/ks以及pn/ps
## kaks需要提取同源的蛋白质以及CDS序列【de novo基因已经有】；需要检查人中的蛋白质序列
## 是否和*ortholog.pep.fa中的匹配；以及需要把首个氨基酸替换成M再进行匹配的比较
## pnps可以直接根据基因组的坐标计算

##### pnps #####
# 1.1 注释vcf
nohup bash evolution_para/1.annotate_vcf.sh > ../log/annotate_vcf.log 2>&1 &
# nohup bash evolution_para/1.annotate_vcf.test.sh > ../log/annotate_vcf.test.log 2>&1 &
# tail -f ../log/annotate_vcf.log

##### kaks #####
# 1.1 单个case测试
## 需要修改paraAT的脚本
denovo_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251029_denovo_check/processed
orf=PB.29415.11__chr22__-__201__3603__2541__2736__dORF__CTG
bash ./kaks/kaks.bin.sh \
    -p $denovo_dir/check_denovo/chr22/orfs/$orf.ortholog.pep.fa \
    -c $denovo_dir/check_denovo/chr22/orfs/$orf.ortholog.nucl.fa \
    -o ../processed/evolution_para/kaks/$orf \
    -m human \
    -t 10
# 1.2 批量运行和合并
## 1.2.1 小批量试运行
find $denovo_dir -name "*ortholog.pep.fa" | head -n50 | while read pep; do
    nucl="${pep%.pep.fa}.nucl.fa"
    cp "$pep" tmp/
    cp "$nucl" tmp/
done
nohup bash ./kaks/kaks.parallel.sh $PWD/tmp ../processed/evolution_para/kaks &> ../log/kaks.parallel.log &
## 1.2.2 正式运行
nohup bash ./kaks/kaks.parallel.sh $denovo_dir ../processed/evolution_para/kaks &> ../log/kaks.parallel.log &
# 监控：
# ls /home/user/data3/lit/project/sORFs/11-denovo-list/processed/evolution_para/kaks|wc -l