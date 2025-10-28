readRDS("./output_20250407/all_orfs_ORF_id_seq_trans_map_correct.rds") -> all_orfs_ORF_id_seq_trans_map_correct
readRDS("./output_20250407/all_orfs_ORF_id_seq_trans_map_incorrect.rds") -> all_orfs_ORF_id_seq_trans_map_incorrect
class(all_orfs_ORF_id_seq_trans_map_correct)
class(all_orfs_ORF_id_seq_trans_map_incorrect)
merge.data.table(all_orfs_ORF_id_seq_trans_map_correct,
                 all_orfs_ORF_id_seq_trans_map_incorrect,by="ORF_id_seq") -> correct_incorrect_map
# 80%的id都是可以map上的
sum(correct_incorrect_map$ORF_id_trans.x==correct_incorrect_map$ORF_id_trans.y)/nrow(correct_incorrect_map)
colnames(correct_incorrect_map) <- c("ORF_id_seq","ORF_id_trans_correct","ORF_id_trans_incorrect")
saveRDS(correct_incorrect_map,"./output_20250407/correct_incorrect_map.rds")
