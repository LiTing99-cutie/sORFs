source activate hla
mkdir -p optitype_out && cd optitype_out
SAMPLE="L1EJF1602305-p21"
BAM=/home/user/data3/lit/project/sORFs/07-Genome/processed/human_brain_21pcw.sorted.bam
[ -f %{BAM}.bai ] || samtools index -@ 40 $BAM
samtools view -@ 40 -b ${BAM} chr6:28000000-34000000 > ${SAMPLE}_hla.bam
samtools sort -@ 40 -n ${SAMPLE}_hla.bam -o ${SAMPLE}_hla_sorted.bam
samtools fastq -1 ${SAMPLE}_hla_R1.fq -2 ${SAMPLE}_hla_R2.fq ${SAMPLE}_hla_sorted.bam
gzip ${SAMPLE}_hla_R*.fq
mkdir -p optitype_${SAMPLE}
OptiTypePipeline.py \
    -i ${SAMPLE}_hla_R1.fq.gz ${SAMPLE}_hla_R2.fq.gz \
    --dna \
    --outdir optitype_${SAMPLE} \
    --verbose