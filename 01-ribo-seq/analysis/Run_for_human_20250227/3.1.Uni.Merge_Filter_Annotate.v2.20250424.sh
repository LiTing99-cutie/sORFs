# refer to /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.1.Uni.Merge_Filter_Annotate.v1.sh

# organize_res_path=$PWD/human_brain_ribo_merge_call_orfs_20250338/organized
# output_path=$PWD/human_brain_ribo_merge_call_orfs_20250338/merge
organize_res_path=$1
output_path=$2
mkdir -p $output_path && cd $output_path
# PRICE
awk -v OFS='\t' '{print $1,$2,"PRICE"}' $organize_res_path/PRICE/nonCano.sorf.tab > nonCano.sorf.meta.raw.PRICE.txt
# RiboCode
awk -v OFS='\t' '{print $1,$2,"RiboCode"}' $organize_res_path/RiboCode/nonCano.sorf.tab > nonCano.sorf.meta.raw.RiboCode.txt
# RibORF
awk -v OFS='\t' '{print $1,$2,"RibORF"}' $organize_res_path/RibORF/nonCano.sorf.tab > nonCano.sorf.meta.raw.RibORF.txt
# merge
cat nonCano.sorf.meta.raw.PRICE.txt nonCano.sorf.meta.raw.RiboCode.txt nonCano.sorf.meta.raw.RibORF.txt | awk -v OFS='\t' '{print $1,$2,$3}' > nonCano.sorf.meta.merge.raw.3_ways.txt
