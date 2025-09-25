cd output
awk -v OFS='\t' '{print $1,$2+1,$3,$4,$5,$6}' sep.add_anno.cds.bed > sep.add_anno.cds.1-based.bed
define_annotation_gencode_v41_human
bash /home/user/data3/lit/project/sORFs/01-ribo-seq/S3.0c.Uni.translate_gtf.v2.20250325.sh sep.add_anno.gtf $fa
# sep.add_anno.cds.1-based.bed
# prot.fa