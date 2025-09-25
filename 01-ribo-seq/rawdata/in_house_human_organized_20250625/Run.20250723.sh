for file in p21_0626_{1..15};do
file1=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/data-20250714/rawdata/$file.R1.raw.fastq.gz
file2=/home/user/data3/lit/project/sORFs/data-download-20250721/MJ20250701254-MJ-D-20250629001-李春琼-纯文库-6个样本/rawdata/$file.R1.raw.fastq.gz
cat $file1 $file2 > $file.fq.gz
done