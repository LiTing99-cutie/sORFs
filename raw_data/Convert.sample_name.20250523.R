# 将原始的metadata表格进一步整理成统一的格式
setwd("/rd1/user/lit/project/sORFs/raw_data")
openxlsx::read.xlsx("./LiChunqiong-MS-49samples-20250512-to_XiaoQi.xlsx") -> meta_20250523
meta_20250523 -> df
# 首先提取酶解方法部分（最后一个下划线后的内容）
df$enzyme_method <- sub(".*_([^_]+)$", "\\1", df$样本名称)
df$Enzyme <- df$enzyme_method
df$Enrichment <- word(df$样本名称,4,sep="_")
# 定义替换规则
replacements <- c(
  "LC-T" = "Trypsin_LysC",
  "LN-T" = "Trypsin_LysN",
  "T-T" = "Trypsin",
  "T-C" = "Trypsin_Chymotrypsin",
  "T-GC" = "Trypsin_GluC",
  "T-AN" = "Trypsin_AspN",
  "T-AC" = "Trypsin_ArgC"
)

# 使用stringr包的str_replace_all进行替换
library(stringr)
df$Enzyme <- str_replace_all(df$Enzyme, replacements)

df$样本名称 %>% str_replace(.,"BSEP_","") %>% str_replace(.,"-","_") -> df$Sample
df$Period <- "21pcw"
df$Label_num <- str_replace(df$装瓶后瓶子上的标记,"BSEP_","")
df %>% select(Sample,Enzyme,Period,Enrichment,Label_num) -> df_1
#导出
write.table(df_1,"./MS_metadata_20250523_human.txt",quote = F,row.names = F,sep = ",",col.names = F)

