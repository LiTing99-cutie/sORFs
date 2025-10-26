# refer to /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1.Uni.Merge_Filter_Annotate.v1.sh

# organize_res_path=$PWD/human_brain_ribo_merge_call_orfs_20250338/organized
# output_path=$PWD/human_brain_ribo_merge_call_orfs_20250338/merge
organize_res_path=$1
output_path=$2
mkdir -p $output_path
# PRICE
awk -v OFS='\t' '{print $1,$2,"PRICE"}' $organize_res_path/PRICE/nonCano.sorf.tab > $output_path/nonCano.sorf.meta.raw.PRICE.txt
# RiboCode
awk -v OFS='\t' '{print $1,$2,"RiboCode"}' $organize_res_path/RiboCode/nonCano.sorf.tab > $output_path/nonCano.sorf.meta.raw.RiboCode.txt
# RibORF
awk -v OFS='\t' '{print $1,$2,"RibORF"}' $organize_res_path/RibORF/nonCano.sorf.tab > $output_path/nonCano.sorf.meta.raw.RibORF.txt
# merge
cat $output_path/nonCano.sorf.meta.raw.PRICE.txt $output_path/nonCano.sorf.meta.raw.RiboCode.txt $output_path/nonCano.sorf.meta.raw.RibORF.txt | awk -v OFS='\t' '{print $1,$2,$3}' > $output_path/nonCano.sorf.meta.merge.raw.3_ways.txt
awk -F'\t' 'BEGIN{OFS="\t"}{
  if (match($1, /[+-]chr[^[:space:]]*/)) {
    s = substr($1, RSTART, RLENGTH)      # 提取 +chr... 或 -chr...
    print $0, s ":" $2                   # 追加新列：±chr...:PEPTIDE
  } else {
    print $0, ""                          # 若未匹配到，就追加空值
  }
}' $output_path/nonCano.sorf.meta.merge.raw.3_ways.txt > $output_path/orfs.3_ways.txt
