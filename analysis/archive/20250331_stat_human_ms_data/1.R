# 整合所有的元数据，加上ms2 count的数据

output_path <- "/rd1/user/lit/project/sORFs/analysis/20250331_stat_human_ms_data/output"
if(!dir.exists(output_path)){
  dir.create(output_path,recursive = T)
}
# 导入数据
fread("/rd1/user/lit/project/sORFs/raw_data/MS_metadata_20250326_new_batch_human.csv",data.table = F) -> sample_metadata
fread("/rd1/user/lit/project/sORFs/output/MS/stat/sample_ms2_count_formal_hm.add_h.txt",data.table = F) -> sample_ms2_count
# 命名列
sample_metadata <- mutate(sample_metadata,V6=NULL)
colnames(sample_metadata) <- c("Sample","Eyzyme","Period","Enrichment","Page_Class")
# 合并ms2 count
merge(sample_metadata,sample_ms2_count) -> sample_metadata
# 添加replicate 列
sample_metadata$Sample %>% str_split("_") %>% do.call(rbind,.) %>% .[,2] -> replicates
sample_metadata$Replicate <- replicates
# fragpipe的输出把所有的-变成了_，为了样本名的统一，也修改
sample_metadata$Sample %>% gsub("-","_",.) -> sample_metadata$Sample
sample_metadata %>% arrange(Enrichment) -> sample_metadata_ordered
sample_metadata_ordered$Sample -> sample_order
fwrite_c(sample_metadata_ordered[,"Sample",drop=FALSE],o("sample_order.txt"))
fwrite_c(sample_metadata_ordered,o("sample_metadata_ordered.txt"))

