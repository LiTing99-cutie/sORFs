res_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251013_denovo_check/processed/check_denovo
mkdir -p "$res_path/merged"

# 合并 ancestors（只保留第一个文件的表头）
find "$res_path" -path "$res_path/merged" -prune -o -type f -name "ancestors" -print0 \
| sort -z \
| xargs -0 awk 'FNR==1 && NR!=1 {next} {print}' \
> "$res_path/merged/merged_ancestors"

# 合并 spec.out（只保留第一个文件的表头）
find "$res_path" -path "$res_path/merged" -prune -o -type f -name "spec.out" -print0 \
| sort -z \
| xargs -0 awk 'FNR==1 && NR!=1 {next} {print}' \
> "$res_path/merged/merged_spec.out"

map_path=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251013_denovo_check/processed/work_mafsInRegion
find $map_path -name "*_safe_map.tsv"|xargs cat > $map_path/safe_map_merged.tsv
