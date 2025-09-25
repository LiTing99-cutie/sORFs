# iso-seq results
iso_res=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/processed/classify/collapsed_classification.filtered_lite_classification.txt
less $iso_res |cut -f 1,6,7,15,22,23 > ../processed/feature_preprare/isoform.expr.info.txt
less $iso_res |cut -f 1,22,23 > ../processed/feature_preprare/isoform.expr.txt
# ribo-seq results
ribo_res=../processed/ribo/counts.txt
# # grep -v ^# $ribo_res|cut -f 1,7,8|awk 'NR==1{print "ORF_id\tRPF_reads\tPsites_number"; next}1' > ../processed/feature_preprare/orf.rpf.psite.txt
# # grep -v ^# $ribo_res|tail -n +2|awk '{print $1,$7/$6*100,$8/$6*100}'| \
# #     awk 'NR==1{print "ORF_id\tRPF_reads\tPsites_number"}1' > ../processed/feature_preprare/orf.rpf.psite.txt
# out=../processed/feature_preprare/orf.rpf.psite.txt
# # 计算总reads（用于“每百万”）
# read rpf_total ps_total < <(
#   grep -v '^#' "$ribo_res" | awk 'NR==1{next} {rpf+=$7; ps+=$8} END{print rpf, ps}'
# )
# # 容错：避免除零
# : "${rpf_total:=0}"
# : "${ps_total:=0}"
# grep -v '^#' "$ribo_res" | awk -v rpfT="$rpf_total" -v psT="$ps_total" '
# BEGIN{
#   OFS="\t";
#   print "ORF_id","RPF_reads","Psites_number","RPF_RPKM","Psites_RPKM"
# }
# NR==1{next}  # 跳过表头
# {
#   id=$1; len_bp=$6+0; rpf=$7+0; ps=$8+0;

#   # RPKM = count / (length_kb * total_million) = count * 1e9 / (len_bp * total_counts)
#   rpf_rpkm = (len_bp>0 && rpfT>0) ? (rpf*1e9)/(len_bp*rpfT) : "NA";
#   ps_rpkm  = (len_bp>0 && psT>0)  ? (ps *1e9)/(len_bp*psT)  : "NA";

#   print id, rpf, ps, rpf_rpkm, ps_rpkm
# }' > "$out"


# ---
out=../processed/feature_preprare/orf.rpf.psite.txt

# 计算总reads（用于“每百万”）
read rpf_total ps_total < <(
  grep -v '^#' "$ribo_res" | awk 'NR==1{next} {rpf+=$7; ps+=$8} END{print rpf, ps}'
)
: "${rpf_total:=0}"
: "${ps_total:=0}"

grep -v '^#' "$ribo_res" | awk -v rpfT="$rpf_total" -v psT="$ps_total" '
BEGIN{
  OFS="\t";
  print "ORF_id","RPF_reads","Psites_number","RPF_RPKM","Psites_RPKM","RPF_codon_coverage","Psites_codon_coverage"
}
NR==1{next}
{
  id=$1; len_bp=$6+0; rpf=$7+0; ps=$8+0;

  # RPKM
  rpf_rpkm = (len_bp>0 && rpfT>0) ? (rpf*1e9)/(len_bp*rpfT) : "NA";
  ps_rpkm  = (len_bp>0 && psT>0)  ? (ps *1e9)/(len_bp*psT)  : "NA";

  # 每密码子覆盖（len_bp 为 bp，1 个密码子 = 3 bp）
  rpf_cov = (len_bp>0) ? (rpf*3)/len_bp : "NA";
  ps_cov  = (len_bp>0) ? (ps *3)/len_bp : "NA";

  print id, rpf, ps, rpf_rpkm, ps_rpkm, rpf_cov, ps_cov
}' > "$out"
