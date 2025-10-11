conda activate base

# 按照样本拆分输入
python3 split_by_sample.py \
  --pep_orf_merged ../MS_res_from_Galaxy/pep.orf.merged.txt \
  --intensity_merged ../MS_res_from_Galaxy/peptide_intensity_IL.merged.tsv \
  --out_base ../processed/by_sample \
  --sample_list ../processed/by_sample/sample_list.txt

# 按照样本运行脚本
BASE=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/log
nohup cat ../processed/by_sample/sample_list.txt| xargs -P 8 -I{} bash -lc 'd="'"$BASE"'/{}/logs"; mkdir -p "$d";bash ./run.all.20250917.v2.sh "{}">"$d/driver.1.log" 2>&1' &

# 合并ibaq_b_with_total结果
python3 merge_add_sample.py \
  --root "../processed/by_sample" \
  --name "ibaq_b_with_total.tsv" \
  --out "../results/all_samples.ibaq_b_with_total.tsv" \
  --mode dirnameN \
  --depth 2

# 折叠ORF，并计算ORF被支持次数
python3 orf_fold_and_occurrence.py \
  --in ../results/all_samples.ibaq_b_with_total.tsv \
  --out ../results/ibaq_orf_folded.tsv

# 合并assignments.unique.post结果
python3 merge_add_sample.py \
  --root "../processed/by_sample" \
  --name "assignments.unique.post.tsv" \
  --out "../results/all_samples.assignments.unique.post.tsv" \
  --mode dirnameN \
  --depth 2

# 计算ORF的unique peptide number
in=../results/all_samples.assignments.unique.post.tsv
out=../results/unique_peptide_counts.tsv
awk -F'\t' '
  NR==1{
    hasHeader = (tolower($1) ~ /^peptide/ || tolower($2) ~ /^assigned/); 
    if (hasHeader) next
  }
  {
    p=$1; id=$2;
    if(p=="" || id=="") next;
    key=p "\t" id;
    if(!(key in seen)){ seen[key]=1; cnt[id]++ }
  }
  END{
    print "ORF_id\tUnique_peptide_n";
    for(id in cnt) print id "\t" cnt[id]
  }
' "$in" | sort -k1,1 > "$out"
echo "[OK] 写出：$out"


if(1==0)
  {
  file1=../results/all_samples.assignments.unique.post.tsv   # 列: Peptide Assigned_protein_id reason Sample
  file2=../MS_res_from_Galaxy/merged.peptide.tsv                            # 列: Peptide_I_L_equal Source Source_number Sample
  out=../results/unique_peptide_counts_with_source.tsv

  awk -F'\t' -v OFS='\t' '
    # ---------- 读 文件2：建立 (Peptide_I_L_equal,Sample) 是否含 msfragger_closed 的映射 ----------
    NR==FNR {
      if (FNR==1) { for(i=1;i<=NF;i++) h2[$i]=i; next }
      k = $h2["Peptide_I_L_equal"] "\t" $h2["Sample"]
      if ($h2["Source"] ~ /msfragger_closed/) closed[k]=1
      next
    }

    # ---------- 读 文件1：统计 ----------
    FNR==1 { for(i=1;i<=NF;i++) h1[$i]=i; next }

    {
      p = $h1["Peptide"]
      id= $h1["Assigned_protein_id"]   # ORF_id
      s = $h1["Sample"]
      if (p=="" || id=="") next

      pair = p "\t" id
      # 总 unique（跨样本去重）
      if (!(pair in seen)) { seen[pair]=1; cnt[id]++ }

      # 若该 (Peptide, Sample) 在文件2中标注过 msfragger_closed，则记为 closed-unique
      k = p "\t" s
      if (k in closed && !(pair in seen_closed)) { seen_closed[pair]=1; cnt_closed[id]++ }
    }

    END {
      print "ORF_id","Unique_peptide_n","Unique_peptide_n_msfragger_closed"
      for (id in cnt) {
        c2 = (id in cnt_closed) ? cnt_closed[id] : 0
        print id, cnt[id], c2
      }
    }
  ' "$file2" "$file1" | sort -k1,1 > "$out"

  echo "[OK] 写出：$out"
  }

# 计算ORF的相对iBAQ
# way 1 使用run中total intensity归一化
python3 orf_mean_relative_ibaq.py \
  --in ../results/all_samples.ibaq_b_with_total.tsv \
  --metric iBAQ_A \
  --out ../results/ibaq_orf_mean_relative.tsv
# way 2 使用所有run中都有丰度的蛋白质做归一化
# python3 ribaq_anchor_norm.py \
#   --in ../results/all_samples.ibaq_b_with_total.tsv \
#   --out ../results/ibaq_orf_mean_relative.tsv \
#   --out_anchors ../results/anchors.allruns.tsv \
#   --out_factors ../results/anchor_scale_factors.tsv \
#   --present_by ibaq_pos \
#   --min_anchor_prop 1.0 \
#   --agg mean
#   # --anchor_only_canonical \

# 将ibaq_orf_folded.tsv与ibaq_orf_mean_relative.tsv以及unique_peptide_counts.tsv合并

# 合并所有信息
rpf=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt
iso=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt
rna=/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt
gene_anno=/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt
python3 merge_tables_all_sample.py \
  --ibaq ../results/ibaq_orf_folded.tsv \
  --rpf  $rpf \
  --iso  $iso \
  --rna  $rna \
  --gene_anno $gene_anno \
  --uniq_pep ../results/unique_peptide_counts_with_source.tsv \
  --rel_ibaq ../results/ibaq_orf_mean_relative.tsv \
  --out ../results/orfs_merged_final.tsv


python3 spearman_corr.v2.py \
  --merged ../results/orfs_merged_final.tsv \
  --out ../results/spearman_results.tsv \
  --pairs A:RPF_RPKM A:mean_relative_iBAQ RPF_RPKM:mean_relative_iBAQ
