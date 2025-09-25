# 2、组装转录本
ref_gtf=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
cd human_brain_rna_seq_alignment_assemble
samtools merge -@ 30 -o merged.bam ./*Aligned.sortedByCoord.out.bam
samtools sort -@ 30 merged.bam -o merged.sorted.bam
stringtie merged.sorted.bam \
          -o stringtie_output/assembled.gtf \
          -p 30 \
          -G $ref_gtf