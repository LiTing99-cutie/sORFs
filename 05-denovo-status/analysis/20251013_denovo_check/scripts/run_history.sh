### 测试 ###

# python3 extract_align_parallel.py --per-chr-bed-root ../processed/get_input/per_chr/ORFs_bed \
#     --maf-root /home/user/data3/lit/project/sORFs/05-denovo-status/maf \
#     --output-root ../processed/runs \
#     --chr chr1 \
#     --workers 5 \
#     --focal hg38

python3 extract_with_mafindex.py \
  --per-chr-bed-root ../processed/get_input/per_chr/ORFs_bed \
  --maf-root /home/user/data3/lit/project/sORFs/05-denovo-status/maf \
  --output-root ../processed/runs_mafindex \
  --work-dir   ../processed/work_mafindex \
  --chr chr1 \
  --focal hg38 \
  --workers 5

# run_mafsInRegion_planA.1.py更名为1_extract_ma_with_mafsInRegion.py
python3 run_mafsInRegion_planA.1.py \
  --per-chr-bed-dir ../processed/get_input/per_chr/ORFs_bed/chr1 \
  --maf-root /home/user/data3/lit/project/sORFs/05-denovo-status/maf \
  --output-root ../processed/runs_mafsInRegion \
  --work-dir   ../processed/work_mafsInRegion \
  --focal hg38

# 重复之前的流程，看结果是否一致【结果是一致的】
scriptDir=/home/user/data3/lit/project/sORFs/05-denovo-status/Denovo_genes-tumors/evolution_orfs
orf_bed_file="/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251013_denovo_check/processed/get_input/per_chr/ORFs_bed/chr1/PB.1.16:chr1:-|66|1611:871:1024|noncoding|GTG.ORF.bed"
mkdir tmp
## 同时更改第四列的内容
mv $orf_bed_file tmp/test.bed
maf_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251013_denovo_check/processed/work_mafindex/_maf_unzip
python3 $scriptDir/1_extract_multiple_alignments.py \
    -b tmp/test.bed -m $maf_dir -o tmp -f yes

# 测试，传递参数染色体
nohup bash run.step_2_3.sh chr1 &> ../log/run.step_2_3.chr1.log &

python3 1_extract_ma_with_mafsInRegion.v1.py \
  --per-chr-bed-dir ../processed/get_input/per_chr/ORFs_bed/$chr \
  --maf-root /home/user/data3/lit/project/sORFs/05-denovo-status/maf \
  --output-root ../processed/test_new_ma/check_denovo \
  --work-dir   ../processed/test_new_ma/work_mafsInRegion \
  --focal hg38

### 正式流程 ###
nohup bash run_extract_ma_with_mafsInRegion.sh &> ../log/run_extract_ma_with_mafsInRegion.log &

nohup bash run.step_2_3.run.sh &> ../log/run.step_2_3.run.log &

# chrX 貌似失败了，考虑重新运行
## 修改脚本设置set -eou pipefail
nohup bash run.step_2_3.sh chrX &> ../log/run.step_2_3.chrX.log &
# nohup bash run.step_2_3.sh chrY &> ../log/run.step_2_3.chrY.log &

nohup bash Step_4.sh &> ../log/Step_4.log &
nohup bash Step_4_supple.sh &> ../log/Step_4_supple.log &