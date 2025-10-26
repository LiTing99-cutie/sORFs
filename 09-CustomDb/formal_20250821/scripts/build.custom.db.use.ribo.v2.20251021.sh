source activate base
candidateORF=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa
# 从中筛选出frame0_fraction大于等于0.5 & Psites_codon_coverage大于等于0.1的ORF；此外，去掉那些和经典CDS重叠的非经典ORF
frame_stats=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/peri/psite_frame_stats.v2.tsv
rpf_psite=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt
overlap=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf_overlap_inframe.txt

output_path=../processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021
mkdir -p $output_path
python3 filter_orf_fasta.v2.py \
  --candidate-fa /home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/candidateORF.6aa.long.M.pep.fa \
  --psite-stats /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/peri/psite_frame_stats.v2.tsv \
  --rpf-psite  /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf.rpf.psite.txt \
  --overlap    /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/orf_overlap_inframe.txt \
  --out-fa     $output_path/filtered_candidates.fa \
  --f0-thr 0.5 --cov-thr 0.1 &> ../log/filter_orf_fasta.v2.log
# 测试同一个任务claude Sonnet 4.5的性能【似乎速度要慢一些】
# python filter_orf_fasta.v2.claude.py \
#     --frame_stats $frame_stats \
#     --rpf_psite $rpf_psite \
#     --overlap $overlap \
#     --input_fasta $candidateORF \
#     --output_fasta $output_path/filtered_candidateORF.fa \
#     --output_table $output_path/filtered_ORF_table.tsv

seqkit rmdup -s -t protein -w 0 -j 20 \
    -D $output_path/candidateORF.filtered.M.dup.txt \
    -o $output_path/candidateORF.filtered.rmdup.pep.fa \
    $output_path/filtered_candidates.fa

python pick_representative.20250910.py \
  -i $output_path/candidateORF.filtered.M.dup.txt \
  -e /home/user/data3/lit/project/sORFs/10-feature-egi/processed/feature_preprare/isoform.expr.txt \
  -o $output_path/representative.tsv \
  --seed 42

python3 replace_ids.py \
  --map $output_path/representative.tsv \
  --fasta-in ${output_path}/candidateORF.filtered.rmdup.pep.fa \
  --fasta-out ${output_path}/candidateORF.filtered.rmdup.renamed.pep.fa

echo "--- Count ORF numbers per ORF type --- "
awk '/^>/{ 
    h = substr($0, 2);            # 去掉 >
    n = split(h, a, "\\|");      # 以 | 分割
    t = (n == 5 ? a[n-1] : "Uniprot Reviewed");
    counts[t]++
}
END {
    for (k in counts) print k "\t" counts[k]
}' ${output_path}/candidateORF.filtered.rmdup.renamed.pep.fa | sort -k2,2nr \
    > ${output_path}/candidateORF.filtered.rmdup.renamed.pep.orf_type.txt

contam_fasta=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/contaminant_fasta/2022_JPR_contam.fasta
cat ${output_path}/candidateORF.filtered.rmdup.renamed.pep.fa $contam_fasta > $output_path/candidateORF.filtered.rmdup.renamed.addContam.pep.fa