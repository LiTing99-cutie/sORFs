nohup bash run.1.sh &> ../log/run.1.log &
nohup bash run.1.rerun.sh &> ../log/run.1.rerun.log &
# 修改部分错误脚本内容，重新运行
nohup bash run.1.rerun.sh &> ../log/run.2.rerun.log &
# 【以最后一版为准，存在修改】
mkdir -p ../processed/annotation/RiboCode_annot
cp /home/user/data/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed/processed/annotation/RiboCode_annot ../processed/annotation/
nohup bash run.1.rerun.sh &> ../log/run.3.rerun.log &
nohup bash run.1.rerun.sh &> ../log/run.4.rerun.log &

# 提取P sites
## 241550834
nohup bash gen.p.sites.bam.sh ../processed/filtered_bam/bam.lst \
    ../processed/qual_assess/ribotish/offset.tab.all.txt ../processed/merged_bam/p_sites_bam 40 &> ../log/gen.p.sites.bam.log &
