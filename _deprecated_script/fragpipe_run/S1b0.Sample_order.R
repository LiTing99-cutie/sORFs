# 导入样本元数据
fread("./raw_data/MS_metadata_20241125.csv",data.table = F) -> sample_metadata
head(sample_metadata)
## fragpipe的输出把所有的-变成了_，为了样本名的统一，也修改
sample_metadata$Experiment %>% gsub("-","_",.) -> sample_metadata$Experiment
sample_metadata %>% arrange(Enrichment) -> sample_metadata_ordered
## RecoveryORnot 回收buffer中的蛋白质太少了，因此不统计了
sample_metadata_ordered %>% filter(RecoveryORnot!="Yes") -> sample_metadata_ordered_subset

sample_metadata_ordered_subset$Experiment -> sample_order
sample_order[c(1,4,2,3,5:24)] -> sample_order
saveRDS(sample_order,file = "./output/MS/visual/sample_order.rds")
saveRDS(sample_metadata_ordered_subset,file = "./output/MS/visual/sample_metadata_ordered_subset.rds")
