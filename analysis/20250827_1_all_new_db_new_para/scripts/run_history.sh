nohup bash fragpipe.20250827.sh &> ../log/fragpipe.20250827.log &
# /rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt

### C8 测试 ###
python merge.db.search.res.20250909.py \
  --msf_closed $(find ../processed/db_search_closed -name "psm.tsv"|grep $sample) \
  --msf_open $(find ../processed/db_search_open -name "psm.tsv"|grep $sample) \
  --pfind_closed $(ls ../processed/pFind_res_20250829/closed/*.spectra |grep $sample) \
  --pfind_open $(ls ../processed/pFind_res_20250829/open/*.spectra |grep $sample) \
  --out ../processed/db_search_merge/$sample.tsv

python fold_by_peptide_il.20250909.py -i ../processed/db_search_merge/$sample.tsv -o ../processed/db_search_merge/$sample.peptide.tsv

protein.map.20250909.c8.test.sh

# c8,得到肽段-intensity-source表格
sample=21pcw_1_C8_T_T
python pep_intensity_merge_il.py \
  --closed $(find ../processed/db_search_closed -name "peptide.tsv"|grep $sample) \
  --open $(find ../processed/db_search_open -name "peptide.tsv"|grep $sample) \
  --sample $sample \
  --out ../processed/quant/$sample/peptide_intensity_IL.tsv \
  --agg max

### 所有的样本 ###
# 步骤一：合并db_search的结果；按照Peptide_I_L_equal聚合并并汇总Source；将肽段比对到蛋白质上
bash separate.pfind.res.20250911.uni.sh

# nohup bash run.all.experi.20250912.sh &> ../log/run.all.experi.20250912.log &
## 运行版本
nohup bash run.all.experi.20250912.parallel.sh &> ../log/run.all.experi.20250912.parallel.log &

# 步骤二：所有样本得到肽段-intensity-source表格
nohup bash quant.all.experi.20250912.sh &> ../log/rquant.all.experi.20250912.log &

# 整合所有结果到一个大表格中

## 1、整合肽段来源表格
meta_file=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
out=../processed/db_search_merge/merged.peptide.tsv
dir=../processed/db_search_merge
mkdir -p "$dir"
: > "$out"   # 清空输出

# 逐个样本合并，给每行加上 Sample 列
while IFS= read -r sample_raw; do
  sample=${sample_raw//-/_}
  in="$dir/$sample.peptide.tsv"
  [[ -s "$in" ]] || { echo "[WARN] 缺失或空文件: $in"; continue; }

  if [[ ! -s "$out" ]]; then
    # 第一个文件：保留表头并新增列名 Sample
    awk -v s="$sample" 'BEGIN{OFS="\t"}
      FNR==1 {print $0, "Sample"; next}
             {print $0, s}
    ' "$in" > "$out"
  else
    # 后续文件：跳过各自表头，只追加数据
    awk -v s="$sample" 'BEGIN{OFS="\t"}
      FNR==1 {next}
             {print $0, s}
    ' "$in" >> "$out"
  fi
done < <(cut -d',' -f1 "$meta_file" | sed 's/-/_/g')

echo "[OK] 合并完成：$out"
## 2、整合肽段ORF对应表格
out=../processed/protein_map/pep.orf.merged.txt
printf "Peptide\tProtein\tSample\n" > "$out"
for f in ../processed/protein_map/*/pep.orf.txt; do
sample=$(basename "$(dirname "$f")")
awk -v s="$sample" 'BEGIN{OFS="\t"} {print $1, $2, s}' "$f" >> "$out"
done
## 3、整合肽段-intensity表格
awk 'FNR==1 && NR!=1 {next} 1' ../processed/quant/*/peptide_intensity_IL.tsv \
  > ../processed/quant/peptide_intensity_IL.merged.tsv
