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
# path=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/ortholog_noncano_sample
# python3 kaks/check_codon.py -d $path -p "*.nucl.fa" -r --only-problems
# path=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare/input_for_ortholog_extraction/ortholog_cano_lnc_intergenic_sample
# python3 kaks/check_codon.py -d $path -p "*.nucl.fa" -r --only-problems

nohup bash kaks/cp.noncano.ortholog.res.sh &> ../log/cp.noncano.ortholog.res.log &
nohup bash kaks/cp.ortholog.res.translate.sh &> ../log/cp.ortholog.res.translate.log &
nohup bash kaks/cp.denovo.ortholog.res.sh &> ../log/cp.ortholog.res.translate.log &
# nohup bash kaks/cp.noncano.orfs.r1_m_minus.ortholog.res.sh &> ../log/cp.ortholog.res.translate.log &
nohup bash kaks/kaks.run.2025114.sh &> ../log/kaks.run.20251114.log &
nohup bash kaks/kaks.run.20251128.sh &> ../log/kaks.run.20251128.log &
# 监控：
# ls /home/user/data3/lit/project/sORFs/11-denovo-list/processed/evolution_para/kaks|wc -l