rm -rf /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/merged.bam
cp /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/*Aligned.toTranscriptome.out.bam \
    /home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/ && rm -rf \
    /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_rna_seq_alignment_assemble/*Aligned.toTranscriptome.out.bam 

src=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227
dst=/home/user/data/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_output_20250227
mkdir -p $dst
nohup cp -r $src $dst && rm -rf $src &