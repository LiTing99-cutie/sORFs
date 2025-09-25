cd tRNA
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir tRNA_STARindex --genomeFastaFiles ./hg38.tRNA.fa --genomeSAindexNbases 6
cd ../rRNA/
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir rRNA_STARindex --genomeFastaFiles ./hg38.rRNA.fa --genomeSAindexNbases 8
cd ../snoRNA/
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir snoRNA_STARindex --genomeFastaFiles ./hg38.snoRNA.fa --genomeSAindexNbases 7