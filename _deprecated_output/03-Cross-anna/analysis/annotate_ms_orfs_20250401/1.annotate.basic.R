```{r}
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")
```

```{r}
# 1.加载R包
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
# the file need to be annotated，需要含有ORF_id_trans列
sep_path <- "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/sep.txt"
# output_path
output_path <- "output"
create_path(output_path)
trans_based_sorfs_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/trans_based_sorfs.txt"
trans_mul_feature_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/Ensembl_106_Gencode_v41_Human_Transcript_stable_ID_version_Gene_stable_ID_version_Gene_name_Transcript_type_gene_type.txt"
gpe_path <- "/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.10.gpe"
correct_incorrect_map_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/correct_incorrect_map.rds"

# 2.整合染色体，起始位置，终止位置，链，长度等信息
fread(sep_path) -> all_orfs
## 修正下ORF_trans_id
correct_incorrect_map <- readRDS(correct_incorrect_map_path)
merge(correct_incorrect_map,all_orfs,by.x="ORF_id_trans_incorrect",by.y="ORF_id_trans") -> all_orfs_new_orf_id_trans
all_orfs_new_orf_id_trans %>% mutate(ORF_id_trans_incorrect=NULL) %>% rename(ORF_id_trans=ORF_id_trans_correct) -> all_orfs
matches <- str_match(all_orfs$ORF_id_trans, "(ENS\\w+\\.\\d+)([+-])(chr\\w+):(\\d+)-(\\d+)")
all_orfs$Strand <- matches[,3]
all_orfs$Chr <- matches[,4]
all_orfs$Start <- matches[,5] %>% as.numeric()
all_orfs$End <- matches[,6] %>% as.numeric()
trans_based_sorfs <- fread_c(trans_based_sorfs_path)
setDT(all_orfs)
setDT(trans_based_sorfs)
setkey(all_orfs, ORF_id_trans)
setkey(trans_based_sorfs, ORF_id_trans)
all_orfs_1 <- trans_based_sorfs[, .(ORF_id_trans, ORF_id_seq, Seq,ENS_id, Scodon, Type)][all_orfs]
all_orfs_1$Length <- nchar(all_orfs_1$Seq)

## 3.整合基因名和转录本类型
trans_mul_feature <- fread(trans_mul_feature_path)
colnames(trans_mul_feature) <- c("Gene_ID","Transcript_ID","Transcript_type","Gene_type","Gene_name")
merge(all_orfs_1,trans_mul_feature[,c("Transcript_ID","Gene_name","Transcript_type","Gene_type")],by.x = "ENS_id",by.y="Transcript_ID") -> all_orfs_1
all_orfs_1 %>% as.data.frame() -> all_orfs_1

# 4.注释小肽类型
gpe <- fread(gpe_path,header = F,data.table = F)
gpe[,c("V1","V6","V7")] -> gpe_selected
colnames(gpe_selected) <- c("ENS_id","CDS_start","CDS_end")
merge(all_orfs_1,gpe_selected) -> all_orfs_2
get_orf_type_2(all_orfs_2) -> all_orfs_3
fwrite_c(all_orfs_3,o("sep.add_anno.txt"))
```



