setwd("/rd1/user/lit/project/sORFs/analysis/20250404_tmp_from_nanjing")
source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
lib_text()

fread_c("../20250331_stat_human_ms_data/output/total_psm.txt") -> total_psm
fread_c("./nonCano.sorf.filtered.add.orf_id_trans.txt") -> ribo_sorfs
fread_c("./output/merged_rna_counts.txt") -> merged_rna_counts
fread_c("ts.meta.txt") -> ts_meta
fread_c("sep_add_overlap_pro_filter.txt") -> in_house_res
# 1. 首先展开所有的肽段
total_psm %>% filter(grepl("ENS",Protein)) %>% filter(Is.Unique=="FALSE") %>% 
  filter(!grepl("sp",Mapped.Proteins)) -> all_sample_sep_group_psm
all_sample_sep_group_psm$Group_id <- 1:nrow(all_sample_sep_group_psm)
all_sample_sep_group_psm %>% separate_rows(`Mapped.Proteins`, convert = TRUE,sep=", ") -> tmp
tmp %>% mutate(Protein=`Mapped.Proteins`,
               `Protein.ID`=`Mapped.Proteins`,
               Gene=`Mapped.Proteins`,
               `Entry.Name`=`Mapped.Proteins`,
               `Protein.Description`=`Mapped.Proteins`) -> tmp_1
rbind(all_sample_sep_group_psm,tmp_1) -> all_sample_sep_group_psm_unfold
# 2. 根据基因表达与否来确定肽段的归属 
all_sample_sep_group_psm_unfold$ENS_id <- str_extract(all_sample_sep_group_psm_unfold$Protein, "ENST[0-9]+\\.[0-9]+")
merge(merged_rna_counts,ts_meta,by = "Geneid",by.y = "V2") %>% select(V1,Mean_rpkm,Mean_counts) -> ens_id_rpkm
colnames(ens_id_rpkm) <- c("ENS_id","Mean_rpkm","Mean_counts")
## 看下之前鉴定得到的小肽是否都是表达大于0的（是的）
merge(ens_id_rpkm,in_house_res) -> in_house_res_add_rpkm
## 比较
merge(all_sample_sep_group_psm_unfold,ens_id_rpkm,by="ENS_id") -> all_sample_sep_group_psm_unfold_add_rpkm
all_sample_sep_group_psm_unfold_add_rpkm %>% filter(Mean_rpkm>0) %>%  group_by(Group_id) %>%
  filter(n() == 1) %>% ungroup()-> psm_retained_group_based_on_transcript
# 基于此救回了50个小蛋白
psm_retained_group_based_on_transcript %>% count(Protein) %>% nrow()
# 3. 根据Ribo-seq来确定肽段的归属 
merge(all_sample_sep_group_psm_unfold,ribo_sorfs,by.x = "Protein",by.y = "ORF_id_trans") -> all_sample_sep_group_psm_unfold_add_ribo_evi
all_sample_sep_group_psm_unfold_add_ribo_evi %>% group_by(Group_id) %>%
  filter(n() == 1) %>% ungroup() -> psm_retained_group_based_on_ribo
# 基于此救回了162个小蛋白
psm_retained_group_based_on_ribo %>% count(Protein) %>% nrow()
# 4.合并结果
colnames(total_psm) -> colnames
rbind(psm_retained_group_based_on_transcript[,colnames],psm_retained_group_based_on_ribo[,colnames]) -> psm_retained
## 有可能基于transcript和基于ribo的结果不一致，但是其实是一致的
intersect(psm_retained_group_based_on_ribo$Group_id,psm_retained_group_based_on_transcript$Group_id) -> id
filter(psm_retained_group_based_on_ribo,Group_id %in% id) %>% .$Protein %>% n_distinct()
all(filter(psm_retained_group_based_on_ribo,Group_id %in% id) %>% .$Protein==
filter(psm_retained_group_based_on_transcript,Group_id %in% id) %>% .$Protein)
## 如果表达大于0和有Ribo-seq证据的冲突，优先选择表达大于0的【但是冲突应该不大】
psm_retained %>% distinct(Spectrum,.keep_all = T) -> psm_retained_1
n_distinct(psm_retained_1$Protein)
fwrite_c(psm_retained_1,"output/psm_retained_based_on_expr_ribo.txt")
## 1276
n_distinct(c(psm_retained_1$Protein,sep_add_overlap_pro$ORF_id_trans))
