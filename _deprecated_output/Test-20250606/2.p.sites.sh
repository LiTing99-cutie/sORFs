# 生成p site的bam文件

cd test_output_20250606

cat all_offset_tab.curated.txt|awk -v OFS='\t' '{print $1,$2+3}' > offset.correction.parameters.txt
bam=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250606/01-output/call-orfs/p21_0523_1/output/alignment/p21_0523_1_Aligned.sortedByCoord.out.bam
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
sam_input=$(echo $bam|sed 's/.bam//').sam
sam_output=p_sites.sam
bam_output=p_sites.bam
samtools view -@ 40 -h $bam > $sam_input
perl $RibORF_path/offsetCorrect.pl -r $sam_input -p offset.correction.parameters.txt -o $sam_output
samtools view -@ 40 -bh $sam_output > $bam_output