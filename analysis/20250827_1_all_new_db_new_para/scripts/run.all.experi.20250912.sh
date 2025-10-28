source activate base
meta_file=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
# for sample in 21pcw_1_3_30K_LC_T;do
for sample in $(cut -f 1 $meta_file -d ','|sed 's/-/_/g'|head -n1);do
    echo -e "***Processing ${sample} at $(date '+%Y-%m-%d %H:%M:%S')"
    # 合并db_search的结果
    python merge.db.search.res.20250909.py \
    --msf_closed $(find ../processed/db_search_closed -name "psm.tsv"|grep $sample) \
    --msf_open $(find ../processed/db_search_open -name "psm.tsv"|grep $sample) \
    --pfind_closed $(ls ../processed/pFind_res_20250829/closed/*.spectra |grep $sample) \
    --pfind_open $(ls ../processed/pFind_res_20250829/open/*.spectra |grep $sample) \
    --out ../processed/db_search_merge/$sample.tsv
    # 按照Peptide_I_L_equal聚合并并汇总Source
    python fold_by_peptide_il.20250909.py -i ../processed/db_search_merge/$sample.tsv -o ../processed/db_search_merge/$sample.peptide.tsv
    # 将肽段比对到蛋白质上
    bash protein.map.20250909.uni.sh $sample
done
