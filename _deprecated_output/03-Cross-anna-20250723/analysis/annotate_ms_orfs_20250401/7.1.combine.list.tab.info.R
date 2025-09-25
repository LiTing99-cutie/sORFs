```{r}
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
uniprot_tab_info_path <- "/home/user/data3/lit/project/sORFs/02-Mass-spec/human/S3/uniprot.human.sep.tab_info.txt"
putative_sep_undetected_sampled_tab_info_path <- "./output/S7/all_putative_sep_undetected_sampled.tab_info.txt"
sep_ms_detected_sampled_tab_info_path <- "./output/sep_add_basic_ms_ribo_info_group_retained.txt"

fread_c(uniprot_tab_info_path) -> uniprot_tab_info
fread_c(putative_sep_undetected_sampled_tab_info_path) -> putative_sep_undetected_sampled_tab_info
fread_c(sep_ms_detected_sampled_tab_info_path) -> sep_ms_detected_sampled_tab_info
mutate(putative_sep_undetected_sampled_tab_info,
       Unique_psm_n=0,
       Unique_peptide_n=0) -> putative_sep_undetected_sampled_tab_info
colnames(uniprot_tab_info) -> x
colnames(putative_sep_undetected_sampled_tab_info) -> y
colnames(sep_ms_detected_sampled_tab_info) -> z
intersect(x,y) %>% intersect(z) -> common_cols

rbind(uniprot_tab_info[,common_cols],
      sep_ms_detected_sampled_tab_info[,common_cols],
      putative_sep_undetected_sampled_tab_info[,common_cols]) -> combined_df

# 生成四个组别
mutate(combined_df,Source=case_when(
  Unique_psm_n!=0 & Type=="Canonical" ~ "Cano_detected",
  Unique_psm_n==0 & Type=="Canonical" ~ "Cano_undetected",
  Unique_psm_n!=0 & Type!="Canonical" ~ "Uncano_detected",
  Unique_psm_n==0 & Type!="Canonical" ~ "Uncano_undetected_sampled",
)) -> combined_df

mutate(combined_df,Source_1=case_when(
  Source=="Cano_detected" | Source=="Cano_undetected"~ "Cano",
  TRUE ~ "Uncano",
)) -> combined_df
```

```{r}
combined_df %>% select(ORF_id_trans) %>% fwrite("./output/S7/orf.id.txt",col.names = T)
```


# 计算疏水性和等电点
```{r}
library(Peptides)
library(Biostrings)
pI(combined_df$Seq) -> combined_df$pI
hydrophobicity(combined_df$Seq, scale = "KyteDoolittle") -> combined_df$Hydrophobicity
```

# 整合表达量
## 整合自产数据的表达量
```{r}
fread_c("/home/user/data3/lit/project/sORFs/06-RNA-seq/02-output/expr/rpkm_N_C_A.txt") -> gene_expr
colnames(gene_expr) <- c("N","C","Gene_id","A")
ts_meta_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/transcript_meta_output/ts.meta.txt"
fread_c(ts_meta_path) -> ts_meta
ts_meta %>% select(2,3) %>% dplyr::rename(Gene_id=V2,Gene_name=V3) %>% 
  distinct(Gene_id,Gene_name) %>% distinct(Gene_name,.keep_all = T)-> ts_meta_1
merge(gene_expr,ts_meta_1) -> gene_expr_add_gene_name
```

```{r}
merge(combined_df,gene_expr_add_gene_name,by.x="Gene_name",by.y="Gene_name") -> combined_df
```
## 整合public数据的表达量
```{r}
fread_c("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S5/merged_rna_counts.txt") -> gene_expr
merge(ts_meta_1,gene_expr,by.x="Gene_id",by.y="Geneid") -> gene_expr_add_gene_name
gene_expr_add_gene_name %>% select(Gene_name,Mean_rpkm) %>% dplyr::rename(Public_mean_rpkm=Mean_rpkm) %>% merge(combined_df,by="Gene_name") -> combined_df
```
## 两个表达量的相关性
```{r}
combined_df %>% filter(Public_mean_rpkm > 0 & A > 0) -> tmp
cor(tmp$Public_mean_rpkm,tmp$A)
```

# 整合不同层次的RPF（多重比对的RPF，特异性比对的RPF，在特定长度的RPF，不同阈值的三碱基周期性的RPF）
## 整理
```{r}
fread_c("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/in_house_phase_I_data_20250520/output/cnt/ribo_counts/RPF.counts.txt") -> ribo_expr
ribo_expr[,c(1,7,8,9,10,11)] -> ribo_expr_1
colnames(ribo_expr_1)
colnames(ribo_expr_1) <- c("ORF_id_trans","RPF_1","RPF_2","RPF_3","RPF_4","RPF_5")
```
## 合并
```{r}
merge(combined_df,ribo_expr_1) -> combined_df
```
## 查看RNA-seq和Ribo-seq的相关性
```{r}
# colnames(combined_df)

cus <- function(df){
c(cor(df$RPF_1/df$Length,df$A),
cor(df$RPF_2/df$Length,df$A),
cor(df$RPF_3/df$Length,df$A),
cor(df$RPF_4/df$Length,df$A),
cor(df$RPF_5/df$Length,df$A)) -> out
return(out)
}
cus(combined_df)
cus(filter(combined_df,Source_1=="Cano"))
cus(filter(combined_df,Source_1=="Uncano"))

cus(filter(combined_df,RPF_1 > 0 & A > 0))
```


# 表达大于0的里面有多少被检测到
```{r}
combined_df %>% count(Source)
filter(combined_df,A>=0.1) %>% count(Source)
filter(combined_df,A>=1) %>% count(Source)
filter(combined_df,A>=10) %>% count(Source)
filter(combined_df,C>=0.1) %>% count(Source)
```
# Ribo-seq大于0的里面有多少被检测到
```{r}
# combined_df %>% count(Source)
cus <- function(cutoff){
print(filter(combined_df,RPF_1>=cutoff) %>% count(Source))
print(filter(combined_df,RPF_2>=cutoff) %>% count(Source))
print(filter(combined_df,RPF_3>=cutoff) %>% count(Source))
print(filter(combined_df,RPF_4>=cutoff) %>% count(Source))
print(filter(combined_df,RPF_5>=cutoff) %>% count(Source))
}
cus(1)
cus(10)
```

# 查看各个类别的数量
```{r}
table(combined_df$Source)
table(combined_df$Source_1)
```
# 查看注释的SEP在多少个基因上
```{r}
filter(combined_df,Source_1=="Cano") -> cano_df
cano_df %>% count(Gene_name) -> cano_df_gene_n
table(cano_df_gene_n$n)
2249/n_distinct(cano_df$Gene_name)
```

# cds_variant被质谱支持的概率更大，其他所有的类型都要更低（需要进一步抽样来判断）
```{r}
combined_df %>% group_by(Source,ORF_type_1) %>% 
  filter(ORF_type_1!="Canonical") %>% 
  summarise(n=n()) -> df
bar_plot_basic_stack(df,x="Source",y="n",fil_col = "ORF_type_1")+
  scale_fill_discrete_divergingx() -> p
p
ggsave(p,filename = "output/S7/plot/uncano_detected_or_not_orf_type.pdf",height = 5,width = 5)
```
# Uncano_detected比Uncano_undetected_sampled长度更长
```{r}
combined_df %>%  ggplot(., aes(x = Length, fill = Source)) +
  geom_density(alpha = 0.7) +
  theme_3()+
  scale_fill_discrete_divergingx()
```
# 起始密码子类型
```{r}
combined_df %>% group_by(Source,Scodon) %>% 
  filter(ORF_type_1!="Canonical") %>% 
  summarise(n=n()) -> df
df$Scodon <- factor(df$Scodon,levels = c("ATG",setdiff(unique(df$Scodon), "ATG")))
bar_plot_basic_stack(df,x="Source",y="n",fil_col = "Scodon")+
  scale_fill_manual(values = brewer.pal(n=n_distinct(df[,"Scodon"]),name="Set3"))
```
# 被Ribo-seq支持的比例
```{r}
combined_df %>% 
  group_by(Source) %>% 
  summarise(
    n=n(),
    ms_support_n=sum(Unique_psm_n!=0),
    ribo_supported_n = sum(Ribo_evidence == 1),
    ribo_supported_ratio = sum(Ribo_evidence == 1) / n(),
    .groups = "drop"  # 取消分组避免后续问题
  )
```
## 被ribo-seq检测到的和被MS检测到的长度是否有区别
```{r}
combined_df %>% filter(Source_1=="Cano") -> tmp
tmp %>% filter(Ribo_evidence==1) -> tmp_ribo
tmp %>% filter(Unique_psm_n>0) -> tmp_MS
mean(tmp_ribo$Length)
mean(tmp_MS$Length)
median(tmp_ribo$Length)
median(tmp_MS$Length)
```


```{r}
select(combined_df,ORF_id_trans) -> df
fwrite(df,"./output/S7/combined_id.txt")
```

# 表达量
```{r}
combined_df %>% box_violin_plot(.,"Source","Gene_mean_rpkm","Source",log10 = T)+
  theme_3(rotate = T)->p
compare(p,list(c("Uncano_detected","Uncano_undetected_sampled"),
               c("Cano_detected","Cano_undetected"))) -> p_1
p_1

combined_df %>% filter(Source!="Uncano_undetected_sampled"& Source!="Cano_undetected") %>% 
ggplot(.,aes_string(x="ORF_type_1",y="Gene_mean_rpkm",fill="ORF_type_1"))+
    geom_violin(trim = F,bw=0.4)+geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+scale_fill_manual(values = brewer.pal(n=n_distinct(combined_df[,"ORF_type_1"]),name="Set3"))+
  theme_3()+
  coord_flip()+
  scale_y_log10(labels = label_number())
```

# 核糖体测序数据检测到的表达是否也会更高
```{r}
combined_df$Ribo_evidence <- as.character(combined_df$Ribo_evidence)
combined_df %>% box_violin_plot(.,"Ribo_evidence","Gene_mean_rpkm","Ribo_evidence",log10 = T)
```

# 疏水性
```{r}
library(scales)
box_violin_plot(data = combined_df,x = "Source",y="Hydrophobicity",fill_col="Source",log10 = F)+
  theme_3(rotate = T) -> p
compare(p,list(c("Uncano_detected","Uncano_undetected_sampled"),
               c("Cano_detected","Cano_undetected"))) -> p_1
p_1+ylim(-4,5)
```
# 等电点
```{r}
box_violin_plot(data = combined_df,x = "Source",y="pI",fill_col="Source",log10 = F)+
theme_3(rotate = T) -> p
compare(p,list(c("Uncano_detected","Uncano_undetected_sampled"),
               c("Cano_detected","Cano_undetected"))) -> p_1
p_1+ylim(3,18)
min(combined_df$pI)
max(combined_df$pI)
```

```{r}
combined_df %>% filter(Source=="Uncano_detected"&ORF_type_1!="Canonical") %>% 
  box_violin_plot(data = .,x = "ORF_type_1",y="pI",fill_col="ORF_type_1",log10 = F)
combined_df %>% filter(Source=="Uncano_detected"&ORF_type_1!="Canonical") %>% 
  box_violin_plot(data = .,x = "ORF_type_1",y="Hydrophobicity",fill_col="ORF_type_1",log10 = F)
```
# 不同类型的疏水性和等电点
## 脑中的非经典，所有组织的经典
```{r}
combined_df %>% filter(Source!="Uncano_undetected_sampled") %>% 
  box_violin_plot(data = .,x = "ORF_type_1",y="pI",fill_col="ORF_type_1",log10 = F)+coord_flip()
combined_df %>% filter(Source!="Uncano_undetected_sampled") %>% 
  box_violin_plot(data = .,x = "ORF_type_1",y="Hydrophobicity",fill_col="ORF_type_1",log10 = F)+coord_flip()

# 调整小提琴图的带宽参数bw，越小越反映局部特征
combined_df %>% filter(Source!="Uncano_undetected_sampled") %>% 
ggplot(.,aes_string(x="ORF_type_1",y="Hydrophobicity",fill="ORF_type_1"))+
    geom_violin(trim = F,bw=0.3)+geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+scale_fill_manual(values = brewer.pal(n=n_distinct(combined_df[,"ORF_type_1"]),name="Set3"))+
  theme_3()+coord_flip()
```
## 脑中的非经典和经典
```{r}
combined_df %>% filter(Source!="Uncano_undetected_sampled"& Source!="Cano_undetected") %>% 
  box_violin_plot(data = .,x = "ORF_type_1",y="pI",fill_col="ORF_type_1",log10 = F)+coord_flip()
# 调整小提琴图的带宽参数bw，越小越反映局部特征
combined_df %>% filter(Source!="Uncano_undetected_sampled"& Source!="Cano_undetected") %>% 
ggplot(.,aes_string(x="ORF_type_1",y="Hydrophobicity",fill="ORF_type_1"))+
    geom_violin(trim = F,bw=0.3)+geom_boxplot(outlier.shape = NA,fill="white",width=0.2)+scale_fill_manual(values = brewer.pal(n=n_distinct(combined_df[,"ORF_type_1"]),name="Set3"))+
  theme_3()+coord_flip()
```

```{r}
fwrite_c(combined_df,"output/S7/combined_df.txt")
fread_c("output/S7/combined_df.txt") -> combined_df
```

