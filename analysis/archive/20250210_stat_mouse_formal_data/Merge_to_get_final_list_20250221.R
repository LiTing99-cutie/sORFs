# 合并不同的质谱run之间的结果
fread_c("output/stat/final_list.txt") -> defaul_ribo
fread_c("output/stat/default_trans/final_list.txt") -> defaul_trans
fread_c("output/stat/open_ribo/final_list.txt") -> open_ribo
fread_c("output/stat/open_trans/final_list.txt") -> open_trans
venn_plot_n(list(defaul_ribo=defaul_ribo$Protein,defaul_trans=defaul_trans$Protein))
venn_plot_n(list(defaul_ribo=defaul_ribo$Protein,open_ribo=open_ribo$Protein))
venn_plot_n(list(defaul_trans=defaul_trans$Protein,open_trans=open_trans$Protein))

cus_1 <- function(df,label){
  df %>% mutate(label=label) -> df
  return(df)
}
rbind(cus_1(defaul_ribo,"defaul_ribo"),
      cus_1(defaul_trans,"defaul_trans"),
      cus_1(open_ribo,"open_ribo"),
      cus_1(open_trans,"open_trans")) -> merged_df
merged_df$label <- factor(merged_df$label,levels = c("defaul_ribo","defaul_trans","open_ribo","open_trans"))
merged_df %>% distinct(Protein,.keep_all = T) -> merged_df_1
table(merged_df_1$label)

# 导出结果到南京上进行进一步的比较分析
# fwrite_c(merged_df_1,path = "output/merged_sorfs.txt")
