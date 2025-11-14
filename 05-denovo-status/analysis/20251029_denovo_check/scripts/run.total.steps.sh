source activate base
orfs_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf
id_list=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/noncano.id.lst
echo -e "***Running get.input at $(date '+%Y-%m-%d %H:%M:%S')"
bash get.input.20251013.sh $orfs_gtf $id_list
echo -e "***Running extract_ma_with_mafsInRegion at $(date '+%Y-%m-%d %H:%M:%S')"
bash run_extract_ma_with_mafsInRegion.sh &> ../log/run_extract_ma_with_mafsInRegion.log
echo -e "***Running Step_2_3 at $(date '+%Y-%m-%d %H:%M:%S')"
bash run.step_2_3.run.sh &> ../log/run.step_2_3.run.log
echo -e "***Running Step_4 at $(date '+%Y-%m-%d %H:%M:%S')"
bash Step_4.sh &> ../log/Step_4.log
