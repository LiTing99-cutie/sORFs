# 设置工作目录
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()

# 读入之前注释好的最终的基于质谱的sep的list
fread_c("./output/sep_add_basic_ms_ribo_info_group_retained.txt") -> ms_sep
ms_sep %>% filter(ORF_type_1=="ncORF") -> ms_sep_on_lncRNA
unique(ms_sep_on_lncRNA$Gene_name) %>% length()

# 如果是ATG起始的密码子，会有些是同一个终止密码子
ms_sep %>% filter(ORF_type_1=="ncORF") %>% filter(Scodon=="ATG") %>% nrow()
ms_sep %>% filter(ORF_type_1=="ncORF") %>% filter(Scodon=="ATG") %>% n_distinct(.$Gene_name)                                     

GO_KEGG_cus <- function(gene_lst){
  library(clusterProfiler)
  library(KEGG.db)
  enrichGO(gene=gene_lst,
           OrgDb="org.Hs.eg.db",
           ont="ALL",
           keyType = "SYMBOL",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05,
           readable = TRUE) -> enrich_GO
  gene_lst_ENTREZID <-
    bitr(gene_lst,
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = "org.Hs.eg.db")
  enrichKEGG(
    gene = gene_lst_ENTREZID$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    use_internal_data = T
  ) -> enrich_KEGG
  add_Desc <- function(enrich_result) {
    frame = toTable(KEGGPATHID2NAME)
    # add description
    merge(frame,
          enrich_result@result,
          by.x = "path_id",
          by.y = "ID") %>% mutate(Description = NULL) %>%
      dplyr::rename(., Description = path_name, ID = path_id) %>% arrange(pvalue, qvalue) -> tmp
    rownames(tmp) <- tmp$ID
    enrich_KEGG_up_add <- enrich_result
    enrich_KEGG_up_add@result <- tmp
    # convert ENTREZID to gene symbol
    setReadable(enrich_KEGG_up_add, 'org.Hs.eg.db', 'ENTREZID') -> enrich_KEGG_up_add_1
    return(enrich_KEGG_up_add_1)
  } 
  add_Desc(enrich_KEGG) -> enrich_KEGG
  return(list(enrich_GO=enrich_GO,enrich_KEGG=enrich_KEGG))
}

gene_lst <- unique(ms_sep_on_lncRNA$Gene_name)
# 所有的sep所在的非编码转录本对应的基因
GO_KEGG_cus(gene_lst) -> ncRNA_sep_GO_KEGG
dotplot(ncRNA_sep_GO_KEGG$enrich_GO)
dotplot(ncRNA_sep_GO_KEGG$enrich_KEGG)

# 所有的sep所在的非编码转录本对应的基因（但是基因是非编码基因）
ms_sep %>% filter(ORF_type_1=="ncORF"&Gene_type!="protein_coding") -> ms_sep_on_lncRNA_non_pc_gene
gene_lst <- unique(ms_sep_on_lncRNA_non_pc_gene$Gene_name)
# NULL
GO_KEGG_cus(gene_lst) -> ncRNA_sep_non_pc_gene_GO_KEGG

library(clusterProfiler)
enrichGO(gene=gene_lst,
         OrgDb="org.Hs.eg.db",
         ont="ALL",
         keyType = "SYMBOL",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05,
         readable = TRUE) -> enrich_GO
dotplot(enrich_GO)


ms_sep_on_lncRNA_non_pc_gene %>% filter(Scodon=="ATG"&Gene_type=="lncRNA")


colnames(ms_sep)
fread_c("./output/S5/merged_rna_counts.txt") -> expr
mean(expr$Mean_rpkm)
median(expr$Mean_rpkm)
# 转录本的信息（用于确定微蛋白对应的转录本在哪个基因上）
ts_meta_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/transcript_meta_output/ts.meta.txt"
fread_c(ts_meta_path) -> ts_meta
distinct(ts_meta,V2,V3) -> gene_id_name_map
merge(expr,gene_id_name_map,by.x="Geneid",by.y="V2") %>% merge(ms_sep,by.x="V3",by.y="Gene_name") -> ms_sep_expr
mean(ms_sep_expr$Mean_rpkm)
median(ms_sep_expr$Mean_rpkm)
filter(ms_sep_expr,Mean_rpkm>=5) %>% nrow()
filter(ms_sep_expr,Mean_rpkm>=10) %>% nrow()

filter(ms_sep_expr,Gene_type!="protein_coding") %>% nrow()
filter(ms_sep_expr,ORF_type_1=="ncORF"&Gene_type!="protein_coding") %>% nrow()
filter(ms_sep_expr,Mean_rpkm>=5 & ORF_type_1=="ncORF"&Gene_type!="protein_coding") %>% nrow()

filter(ms_sep_expr,ORF_type_1=="ncORF"&Gene_type!="protein_coding") %>% density_plot(.,"Mean_rpkm")+xlim(0,5)
filter(ms_sep_expr,ORF_type_1=="ncORF"&Gene_type!="protein_coding") %>% density_plot(.,"Mean_counts")+xlim(0,100)
