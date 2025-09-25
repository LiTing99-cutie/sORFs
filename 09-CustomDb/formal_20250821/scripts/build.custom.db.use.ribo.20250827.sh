source activate base
# 修改
WORK_DIR=../processed/annotation/RibORF_annot
GENOME_FILE=/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/custom_ref.fa
GTF_FILE=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
GPE_FILE=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.15.gpe
RibORF=/home/user/data3/rbase/opt/RibORF/RibORF.2.0
CANDIDATE_ORF_DIR=$WORK_DIR/candidate_ORFs_ribo_filtering
UNIPROT_FA=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.1.fasta
mkdir -p $CANDIDATE_ORF_DIR

# [optional] choose those transctipts with ribo-seq 25-34 nt reads >= 10
FA=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250827_custom_db_filtering/processed/filtering/final.fa
SCR=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/20250827_custom_db_filtering/scripts/custom_db_filtering.20250827.sh
# [ -f "$FA" ] || ( cd "$(dirname "$SCR")" && bash "$(basename "$SCR")" )
pushd "$(dirname "$SCR")" && bash "$(basename "$SCR")" && popd
cp "$FA" "$CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.pep.fa"

# Add uniprot reviewed proteins in the last
echo "--- Add uniprot reviewed proteins in the last --- "
cat $UNIPROT_FA >> $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.pep.fa

# Unique fasta seq of ORF and save duplicated ORFs
echo "--- Unique fasta seq of ORFs and save duplicated ORFs --- "
# -s, --by-seq           by seq
# -w, --line-width int   line width when outputting FASTA format (0 for no wrap) (default 60)
# -o, --out-file string  out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
# -t, --seq-type string  sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
# -j, --threads int      number of CPUs. can also set with environment variable
[ -f $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.rmdup.pep.fa ] || seqkit rmdup -s -t protein -w 0 -j 20 \
    -D $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.dup.txt \
    -o $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.rmdup.pep.fa \
    $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.pep.fa
# first item of candidateORF.6aa.long.M.dup.txt was retained in candidateORF.6aa.long.M.rmdup.pep.fa

# stat3
echo "--- Count ORF numbers per ORF type --- "
awk '/^>/{ 
    h = substr($0, 2);            # 去掉 >
    n = split(h, a, "\\|");      # 以 | 分割
    t = (n == 5 ? a[n-1] : "Uniprot Reviewed");
    counts[t]++
}
END {
    for (k in counts) print k "\t" counts[k]
}' $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.rmdup.pep.fa | sort -k2,2nr \
    > $CANDIDATE_ORF_DIR/candidateORF.6aa.long.M.rmdup.orf_type.txt

# evaluate periodicity
# echo "--- Evaluate periodicity and RRS for all candidate ORFs --- "
# [ -f ./evaluator.log ] || python ./three_nucleotide_periodicity_evaluator_v2.py > ./evaluator.log