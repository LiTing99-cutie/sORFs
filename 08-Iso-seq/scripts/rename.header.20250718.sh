cd /home/user/data3/lit/project/sORFs/08-Iso-seq/rawdata/merge

# bam_1=../r84130_250703_001_1_A01/p21-IsoSeq.Iso_bc02.bcM0004.ISO.bam
input=p21-IsoSeq.Iso_bc02.bcM0004.ISO.merged.20250716.bam
output=p21-IsoSeq.Iso_bc02.bcM0004.ISO.merged.20250718.reNameHeader.bam
samtools view -H $input | sed 's/SM:BioSample_10/SM:BioSample_6/g' |samtools reheader - $input > $output
samtools view -H $output > grep.sample.txt

pbindex -j 30 $output