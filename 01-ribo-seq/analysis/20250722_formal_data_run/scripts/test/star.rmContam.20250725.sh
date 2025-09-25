trimmed_fastq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1/fastqc/p21_0626_2/output/trimmed_fastq/p21_0626_2_trimmed.fq.gz

rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/tmp/blastn/total_rRNA_STARindex
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix ../tmp/rRNA.unalign.20250725.2 \
--outSAMtype None \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1 \
--alignEndsType EndToEnd \
--outFilterMultimapNmax 10000 \
--alignIntronMax 1

rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/tmp/blastn/rRNA_STARindex
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix ../tmp/rRNA.unalign.20250725. \
--outSAMtype None \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1 \
--alignEndsType EndToEnd \
--outFilterMultimapNmax 10000 \
--alignIntronMax 1
