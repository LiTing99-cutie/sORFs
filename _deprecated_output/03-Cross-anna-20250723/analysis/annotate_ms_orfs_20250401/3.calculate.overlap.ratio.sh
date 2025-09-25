source activate base
output_path=/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output
annotation_path=/home/user/data2/lit/project/ZNF271/data/annotation
sep=$output_path/sep.add_anno.txt
sep_gpe=$output_path/sep.add_anno.gpe
sep_gtf=$output_path/sep.add_anno.gtf
sep_cds_bed=$output_path/sep.add_anno.cds.bed
gpe=$annotation_path/gencode.v41.annotation.10.gpe
gtf=$annotation_path/gencode.v41.annotation.gtf
anno_cds_bed=$annotation_path/gencode.v41.annotation.cds.bed
script=/home/user/data3/lit/project/sORFs/03-Cross-anna/Uni.gen.genepred.i_sorf_id.py
python $script $sep $gpe $sep_gpe
genePredToGtf file $sep_gpe $sep_gtf

get_cds_bed(){
	less  $1 | awk '$3=="CDS"'| awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > $2
}

# 得到cds bed
get_cds_bed $sep_gtf $sep_cds_bed
get_cds_bed $gtf $anno_cds_bed
# 求overlap区域
bedtools intersect -a $sep_cds_bed -b $anno_cds_bed -s | sort -k1,1 -k2,2n > $output_path/intersect.sorted.txt
# 每个sep计算重叠区域
cd $output_path
[ -d tmp ] && rm -rf tmp/* 
mkdir -p tmp
for orf_id in $(less intersect.sorted.txt |cut -f 4|sort|uniq);do
	grep $orf_id intersect.sorted.txt > tmp/$orf_id.txt
	bedtools merge -i tmp/$orf_id.txt -c 4,5,6 -o distinct > tmp/$orf_id.dedup.txt
done
cat tmp/*.dedup.txt > intersect.dedup.grouped.txt
awk -v OFS='\t' '{sum[$4] += $3 - $2} END {for (key in sum) print key, sum[key]}' intersect.dedup.grouped.txt > sep.overlapped.nt.txt