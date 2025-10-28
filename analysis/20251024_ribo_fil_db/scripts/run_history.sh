nohup bash fragpipe.20250827.sh &> ../log/fragpipe.20250827.log &
# /rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
# 重新运行的话需要更改pFind文件夹的名字
### 所有的样本 ###
# 步骤一：合并db_search的结果；按照Peptide_I_L_equal聚合并并汇总Source；将肽段比对到蛋白质上
bash separate.pfind.res.20250911.uni.sh

## 运行版本
conda activate base
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
