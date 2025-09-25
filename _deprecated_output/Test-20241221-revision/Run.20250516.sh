script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20241221-revision/Run.raw_reads.qc.20250516.sh
fastq_path=/home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/Test-20241221/

for sample in 110N 110T;do
bash $script $fastq_path/$sample ./output
done

script_3=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250408/Uni.call.orfs.ribocode.human.20250408.sh
find ./ -name "*trimmed.fastq.1.gz"
ln -s ./output/rawdata/110N/output/trimmed_fastq/110N.trimmed.fastq.1.gz ./110N_trimmed.fastq.gz
ln -s ./output/rawdata/110T/output/trimmed_fastq/110T.trimmed.fastq.1.gz ./110T_trimmed.fastq.gz
# call orfs
for trimmed_fq in $(ls $PWD/*_trimmed.fastq.gz);do
bash $script_3 $trimmed_fq $PWD/01-output/call-orfs
done

