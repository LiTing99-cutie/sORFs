# 合并不同的质谱run之间的结果
fread_c("output/stat/final_list.txt") -> defaul_ribo
fread_c("output/stat/default_trans/final_list.txt") -> defaul_trans
fread_c("output/stat/open_ribo/final_list.txt") -> open_ribo
venn_plot_n(list(defaul_ribo=defaul_ribo$Protein,defaul_trans=defaul_trans$Protein))
venn_plot_n(list(defaul_ribo=defaul_ribo$Protein,open_ribo=open_ribo$Protein))

cus_1 <- function(df,label){
  df %>% mutate(label=label) -> df
  return(df)
}
rbind(cus_1(defaul_ribo,"defaul_ribo"),
      cus_1(defaul_trans,"defaul_trans"),
      cus_1(open_ribo,"open_ribo")) -> merged_df
merged_df$label <- factor(merged_df$label,levels = c("defaul_ribo","defaul_trans","open_ribo"))
merged_df %>% distinct(Protein,.keep_all = T) -> merged_df_1
table(merged_df_1$label)

# 转录组鉴定到的有多少是转录组独有的
fread("./sorfs_id/trans_based.sorf.id.txt",data.table = F,header = F) -> trans_based_sorf_id
fread("./sorfs_id/ribo_based.sorf.id.txt",data.table = F,header = F) -> ribo_based_sorf_id
venn_plot_n(list(trans_based_sorf=trans_based_sorf_id$V1,ribo_based_sorf=ribo_based_sorf_id$V1))

# 95/1898 是来自于ribo-seq鉴定到的，转录组独特的有1803个
filter(merged_df_1,label=="defaul_trans")$Protein %in% ribo_based_sorf_id$V1 %>% sum()

# open_ribo鉴定到的有多少是带有修饰的
fread_c("output/stat/open_ribo/all_sample_all_level_sep.txt") -> all_sample_all_level_sep_open_ribo
all_sample_all_level_sep_open_ribo %>% filter(Confidence_Level==1) -> df_level_1
merged_df_1 %>% filter(label=="open_ribo") -> tmp
df_level_1[df_level_1$Protein %in% tmp$Protein,] -> res
res$`Razor Observed Modifications` != "" -> res$modi_or_not
# 从蛋白质的结果来看，只有34/426存在修饰；之前基于一个样本的结果，PSM中有30%多存在修饰，可能没有反映到蛋白质上？
res %>% group_by(Protein) %>% summarise(modi_stat_in_all_sam=sum(modi_or_not)) %>% 
  filter(modi_stat_in_all_sam!=0) %>% nrow()

# 导出结果到南京上进行进一步的比较分析
fwrite_c(merged_df_1,path = "output/merged_sorfs.txt")
