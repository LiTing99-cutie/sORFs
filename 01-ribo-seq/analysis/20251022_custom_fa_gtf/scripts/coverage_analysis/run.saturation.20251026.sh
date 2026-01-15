gc_out_dir=../processed/geneCounts/
res_out_dir=../results/coverage_analysis
mkdir -p $gc_out_dir $res_out_dir
## 计算reads个数
bash coverage_analysis/Uni.saturation.calcu.gene.counts.v1.20250728.sh ../processed/filtered_bam $gc_out_dir &> ../log/saturation.calcu.gene.counts.20250728.log
## 饱和分析
bash coverage_analysis/Uni.transcript_coverage_analysis_v1_20250725.sh $gc_out_dir/transcript.counts.txt  $res_out_dir/transcript
bash coverage_analysis/Uni.transcript_coverage_analysis_v1_20250725.sh $gc_out_dir/gene.counts.txt  $res_out_dir/gene

Rscript coverage_analysis/plot_transcript_coverage_combine_v2_20251026.R $res_out_dir/transcript/transcript_coverage_results.csv \
    $res_out_dir/gene/transcript_coverage_results.csv $res_out_dir/combine.pdf

nohup bash coverage_analysis/run.orf.identification.saturation.20251225.sh &> ../run.orf.identification.saturation.20251225.log &

python coverage_analysis/plot_stats.py