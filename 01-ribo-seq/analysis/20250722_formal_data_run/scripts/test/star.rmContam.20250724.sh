trimmed_fastq=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1/fastqc/p21_0626_2/output/trimmed_fastq/p21_0626_2_trimmed.fq.gz
rRNA_index=/home/user/data3/lit/project/sORFs/01-ribo-seq/Pre-Run/ncRNA/STARindex/rRNA/rRNA_STARindex
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix ../tmp/rRNA.unalign. \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1 \
--alignEndsType EndToEnd \
--alignIntronMax 1
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix ../tmp/rRNA.unalign.1 \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
STAR --runThreadN 50 \
--limitBAMsortRAM 250000000000 \
--seedSearchStartLmax 22 \
--genomeDir $rRNA_index \
--readFilesIn $trimmed_fastq \
--readFilesCommand zcat \
--outFileNamePrefix ../tmp/rRNA.unalign.2 \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--outFilterMismatchNmax 1