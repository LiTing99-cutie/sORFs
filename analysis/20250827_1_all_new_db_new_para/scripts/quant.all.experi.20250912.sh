source activate base
meta_file=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
# for sample in 21pcw_1_3_30K_LC_T;do
for sample in $(cut -f 1 $meta_file -d ','|sed 's/-/_/g');do
    echo -e "***Processing ${sample} at $(date '+%Y-%m-%d %H:%M:%S')"
    python pep_intensity_merge_il.py \
    --closed $(find ../processed/db_search_closed -name "peptide.tsv"|grep $sample) \
    --open $(find ../processed/db_search_open -name "peptide.tsv"|grep $sample) \
    --sample $sample \
    --out ../processed/quant/$sample/peptide_intensity_IL.tsv \
    --agg max
done
