define_annotation_gencode_v41_human
cd output/S9
bedtools getfasta -s -fi $fa -bed window_15_start_codon.bed > window_15_start_codon.fa
for i in window_15_start_codon_cano window_15_start_codon_uncano;do
    bedtools getfasta -s -fi $fa -bed $i.bed > $i.fa
done