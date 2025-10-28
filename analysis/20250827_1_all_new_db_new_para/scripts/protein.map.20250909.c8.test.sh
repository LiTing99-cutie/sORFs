# 第一列 = 已 I→L 的去重肽段
sample=21pcw_1_C8_T_T
out_dir=../processed/protein_map/$sample
mkdir -p $out_dir
cut -f1 ../processed/db_search_merge/$sample.peptide.tsv | tail -n +2 > $out_dir/peptides_il.txt
custom_db=/rd1/user/lit/project/sORFs/custom_database/human/custom_db_20250826_v2/human_brain_custom_db.fasta
pepMatch_jar=/rd1/user/lit/project/sORFs/analysis/archive/20250218_recheck_unique_peptide/PeptideMatchCMD_1.1.jar
java -jar $pepMatch_jar -a index -d $custom_db -i custom_db_index
java -jar $pepMatch_jar -l -a query -i custom_db_index -Q $out_dir/peptides_il.txt -e -o $out_dir/out_fasta.txt &> $out_dir/out.match_n.txt

awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} NF>=5 && length($1)==($5-$4+1)' \
  "$out_dir/out_fasta.txt" > "$out_dir/out_fasta.len_eq.tsv"
  
in="$out_dir/out_fasta.len_eq.tsv" # PeptideMatch结果文件（跳过以#开头的行）
out="$out_dir/peptide_mapped_genes.tsv"

awk -F'\t' '
BEGIN{OFS="\t"}
/^#/ {next}                    # 跳过注释/表头
NF>=2 {
  pep=$1; prot=$2
  if (pep=="" || prot=="") next
  key=pep SUBSEP prot
  # 对每个(肽,蛋白)仅保留一次，且按出现顺序累计
  if (!(key in seen)) {
    seen[key]=1
    if (!(pep in list)) {
      list[pep]=prot
      count[pep]=1
      order[++n]=pep      # 记录肽段首次出现顺序
    } else {
      list[pep]=list[pep]","prot
      count[pep]++
    }
  }
}
END{
  print "Peptide_I_L_equal","Mapped_genes","Mapped_genes_number"
  for (i=1;i<=n;i++){
    p=order[i]
    print p, list[p], count[p]
  }
}' "$in" > "$out"

echo "Wrote: $out"

less "$out_dir/peptide_mapped_genes.tsv" |awk '$3==1'

# 比较和之前方式的区别
## 只是closed search的结果
psm=../processed/db_search_closed/Trypsin/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d/psm.tsv
less ../processed/db_search_closed/Trypsin/CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2_3_1_7020_d/psm.tsv|grep -i true |cut -f 34|grep ^PB|grep -v canonical|sort -u |wc -l
less "$out_dir/peptide_mapped_genes.tsv" |awk '$3==1'|cut -f2|grep ^PB|grep -v canonical|sort -u |wc -l
comm -12 <(less $psm | grep -i true | cut -f34 | grep ^PB | grep -v canonical | sort -u) \
               <(awk -F"\t" '$3==1{print $2}' "$out_dir/peptide_mapped_genes.tsv" | grep ^PB | grep -v canonical | sort -u) \
| wc -l
