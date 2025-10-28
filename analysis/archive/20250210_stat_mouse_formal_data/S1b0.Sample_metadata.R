# 整合所有的元数据，加上ms2 count的数据

# 导入样本元数据
fread("/rd1/user/lit/project/sORFs//raw_data/MS_metadata_20250120_new_batch_mouse.csv",data.table = F) -> sample_metadata
head(sample_metadata)
# 修改元数据中的错误
sample_metadata$V3 <- "E16"
# 命名列
sample_metadata <- mutate(sample_metadata,V6=NULL)
colnames(sample_metadata) <- c("Experiment","Eyzyme","Period","Enrichment","Page_Class")
## fragpipe的输出把所有的-变成了_，为了样本名的统一，也修改
sample_metadata$Experiment %>% gsub("-","_",.) -> sample_metadata$Experiment
sample_metadata %>% arrange(Enrichment) -> sample_metadata_ordered
sample_metadata_ordered$Experiment -> sample_order
saveRDS(sample_order,file = "./output/sample_order.rds")
saveRDS(sample_metadata_ordered,file = "./output/sample_metadata_ordered.rds")

fread("/rd1/user/lit/project/sORFs/output/MS/stat/sample_ms2_count_formal_hm.add_h.txt",data.table = F) -> sample_ms2_count
sample_ms2_count$Sample %>% gsub("-","_",.) -> sample_ms2_count$Sample

merge(sample_metadata_ordered,sample_ms2_count,by.x="Experiment",by.y="Sample") -> sample_metadata
saveRDS(sample_metadata,file = "./output/sample_metadata.rds")

