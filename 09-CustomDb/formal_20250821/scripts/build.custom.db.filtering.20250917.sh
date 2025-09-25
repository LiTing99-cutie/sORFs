source activate base
out_dir=../processed/annotation/RibORF_annot/candidate_ORFs/filtering
UNIPROT_FA=/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/human/uniprotkb_taxonomy_id_9606_AND_reviewed_2025_03_24.1.fasta
mkdir -p $out_dir
python3 filter_orf_fasta.py \
  --rpkm /home/user/data3/lit/project/sORFs/06-RNA-seq/02-output-isoseq-gtf-20250909/expr/rpkm_N_C_A.txt \
  --isoform_map /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.info.txt \
  --ribo /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --fasta /home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.rmdup.pep.fa \
  --out_fasta $out_dir/candidateORF.filtered.fa \
  --trace_tsv $out_dir/candidateORF.filter_trace.tsv \
  --c_thresh 2 \
  --rpf_metric RPF_codon_coverage --ps_metric Psites_codon_coverage \
  --rpf_min 1 --psite_min 1 \
  --evidence_mode any
# Add uniprot reviewed proteins in the last
echo "--- Add uniprot reviewed proteins in the last --- "
cat $UNIPROT_FA >> $out_dir/candidateORF.filtered.fa
# Unique fasta seq of ORF and save duplicated ORFs
echo "--- Unique fasta seq of ORFs and save duplicated ORFs --- "
# -s, --by-seq           by seq
# -w, --line-width int   line width when outputting FASTA format (0 for no wrap) (default 60)
# -o, --out-file string  out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
# -t, --seq-type string  sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
# -j, --threads int      number of CPUs. can also set with environment variable
seqkit rmdup -s -t protein -w 0 -j 20 \
    -D $out_dir/candidateORF.filtered.M.dup.txt \
    -o $out_dir/candidateORF.filtered.rmdup.pep.fa \
    $out_dir/candidateORF.filtered.fa
contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/contaminant_fasta/2022_JPR_contam.fasta
cat $out_dir/candidateORF.filtered.rmdup.pep.fa $contam_fasta > $out_dir/candidateORF.rmdup.pep.addContam.fa