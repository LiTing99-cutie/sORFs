for chr in $(cat ../processed/get_input/chr.list);do
    echo -e "***Processing $chr at $(date '+%Y-%m-%d %H:%M:%S')"
    python3 1_extract_ma_with_mafsInRegion.py \
    --per-chr-bed-dir ../processed/get_input/per_chr/ORFs_bed/$chr \
    --maf-root /home/user/data3/lit/project/sORFs/05-denovo-status/maf \
    --output-root ../processed/check_denovo \
    --work-dir   ../processed/work_mafsInRegion \
    --focal hg38
done