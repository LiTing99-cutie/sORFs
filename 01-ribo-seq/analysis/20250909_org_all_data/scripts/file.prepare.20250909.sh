# 0523_5到0523_7；0626_1到0626_15；0729_1到0729_18
bam_lst_batch_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed/batch_1/filtered_bam/bam.lst
bam_lst_batch_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed/batch_1/filtered_bam/bam.lst
offset_tab_all_1=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/results/batch_1/qual_assess/ribotish/offset.tab.all.txt
offset_tab_all_2=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/results/batch_1/qual_assess/ribotish/offset.tab.all.txt

# 合并
cat $bam_lst_batch_1 $bam_lst_batch_2 > ../results/bam.lst
cat $offset_tab_all_1 $offset_tab_all_2 > ../results/offset.tab.all.txt

less ../results/offset.tab.all.txt|awk -v OFS='\t' '{print $1,$2,$3+3}' > ../processed/offset.tab.all.ribtish.txt

RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
perl $RibORF_path/offsetCorrect.pl -r $sam_input -p offset.correction.parameters.txt -o $sam_output
bam_output=$(echo $sam_output|sed 's/.sam/.bam/')
samtools view -@ 40 -bh $sam_output > $bam_output && rm -rf $sam_output