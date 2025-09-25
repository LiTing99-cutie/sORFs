
# mm39
fa=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
genePred=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.genePred.txt
RibORF_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
perl $RibORF_path/ORFannotate.pl -g $fa -t $genePred -s ATG -l 21 -o annot/RiboORF/mm39
# 增加起始密码子
perl $RibORF_path/ORFannotate.pl -g $fa -t $genePred -s ATG/CTG/ACG/GTG/TTG/ATA/ATC/ATT/AAG/AGG -l 21 -o annot/RiboORF/mm39

# 为blasp建立索引
makeblastdb -in /home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.faa -dbtype prot -out /home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup

# 转换为tab，为了去重
cd /home/user/data3/lit/project/sORFs/01-ribo-seq/annot
# 68205
seqkit seq -s -w0 NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.faa > NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.seq
# 25596
seqkit rmdup -s uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform_2024_10_25.fasta > uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform.rmdup.fasta
seqkit seq -s -w0 uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform.rmdup.fasta > uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform.rmdup.seq
