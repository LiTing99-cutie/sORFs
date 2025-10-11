out=../processed/merged_bam/Aligned.toTranscriptome.out.bam
samtools merge -@ 30 -f -o $out \
  /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed/trans_alignment_change_gtf/*.bam