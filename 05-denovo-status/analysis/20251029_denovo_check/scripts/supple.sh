for chr in $(cat ../processed/get_input/chr.list);do
    echo -e "***Processing $chr at $(date '+%Y-%m-%d %H:%M:%S')"    
    bash run.step_2_3_tmp.sh $chr
done

orfs_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/sub.gtf
id_list=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/noncano.id.lst
total_pep_fa=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/translate_out/prot.fa
total_cds_fa=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20251028_new_db_search_res/results/augment_orf_table/translate_out/cds.fa
nohup bash Step_4.sh $total_pep_fa $total_cds_fa $id_list &> ../log/Step_4.log &