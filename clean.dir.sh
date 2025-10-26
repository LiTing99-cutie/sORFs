rm -rf /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/merged.bam
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/*Aligned.toTranscriptome.out.bam \
    /home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/ && rm -rf \
    /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/*Aligned.toTranscriptome.out.bam 

src=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227
dst=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227
mkdir -p $dst
nohup cp -r $src $dst && rm -rf $src &

# 将小鼠的结果【暂时用不到】转移到data盘中
nohup cp -r /home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011 /home/user/data/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011 &
# 将原始数据转移到data备份盘中
nohup cp -r /home/user/data3/lit/project/sORFs/01-ribo-seq/rawdata/demo_20250813/MJ20250804226-ZX-D-250801018-李春琼-纯文库-18个样本/rawdata/*.R1.raw.fastq.gz /home/user/data/lit/project/sORFs/01-ribo-seq/rawdata/in_house/formal/demo &

# 将比对到非定制化基因组的结果转移到data盘中
src_dir=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed
target_dir=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/20250722_formal_data_run/processed
mkdir -p $target_dir
nohup cp -r $src_dir $target_dir &

#【已经完成】
src_dir=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed
target_dir=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/20250813_demo_data_analysis/processed
mkdir -p $target_dir
nohup cp -r $src_dir $target_dir &