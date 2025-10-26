nohup bash ribo.expr.whole.db.20251024.sh > ../ribo.expr.whole.db.20251024.log &

conda activate base
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20251022_custom_fa_gtf/processed/merged_bam/p_sites_bam/all.offsetCorrected.merged.sorted.bam
gtf=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/ribo/sub.gtf
output_path=../processed/feature_preprare/peri/
mkdir -p $output_path

nohup python3 psite_frame_stats.v2.py \
  --bam $bam \
  --annot $gtf \
  --format gtf \
  --key-attr gene_id \
  --out $output_path/psite_frame_stats.v2.tsv &> ../log/calcu_peri.v2.log &

orf_overlap_inframe_py=orf_overlap_inframe.v2.py
gtf=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/ribo/sub.gtf
anno_gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
mkdir ../processed/feature_preprare/tmp
nohup python3 $orf_overlap_inframe_py \
  --new-gtf $gtf \
  --ann-gtf $anno_gtf \
  --out ../processed/feature_preprare/orf_overlap_inframe.txt \
  --tmpdir ../processed/feature_preprare/tmp \
  --faidx /home/user/data/lit/database/public/genome/hg38/hg38.fa.fai &> ../log/orf_overlap_inframe.v2.log &
