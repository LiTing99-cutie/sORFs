setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human")
source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
# 所有的肽段谱图匹配的路径
total_psm_path <- "./total_psm.txt"
# ribo-seq鉴定出来的微蛋白的路径
ribo_sorfs_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_brain_ribo_merge_call_orfs_20250338/merge/nonCano.sorf.filtered.add.orf_id_trans.txt"
# 微蛋白所在的基因的表达量
rna_seq_counts_path <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/S5/merged_rna_counts.txt"
# 转录本的信息（用于确定微蛋白对应的转录本在哪个基因上）
ts_meta_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/transcript_meta_output/ts.meta.txt"
# 用于查看是否所有MS鉴定出的微蛋白都是表达量大于0的
ms_sorf_path <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/sep_add_overlap_pro.txt"
fread_c(total_psm_path) -> total_psm
fread_c(ribo_sorfs_path) -> ribo_sorfs
fread_c(rna_seq_counts_path) -> merged_rna_counts
fread_c(ts_meta_path) -> ts_meta
fread_c(ms_sorf_path) -> in_house_res
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
## 一共有2207个肽段是比对到多个微蛋白上的
all_sample_sep_group_psm_unfold$Peptide %>% unique() %>% length()
# 2. 根据基因表达与否来确定肽段的归属 
all_sample_sep_group_psm_unfold$ENS_id <- str_extract(all_sample_sep_group_psm_unfold$Protein, "ENST[0-9]+\\.[0-9]+")
merge(merged_rna_counts,ts_meta,by = "Geneid",by.y = "V2") %>% select(V1,Mean_rpkm,Mean_counts) -> ens_id_rpkm
colnames(ens_id_rpkm) <- c("ENS_id","Mean_rpkm","Mean_counts")
## 看下之前鉴定得到的微蛋白是否都是表达大于0的（只有0.3%（3个）的微蛋白表达量为0）
merge(ens_id_rpkm,in_house_res) -> in_house_res_add_rpkm
sum(in_house_res_add_rpkm$Mean_counts==0)/nrow(in_house_res_add_rpkm)
## 比较
merge(all_sample_sep_group_psm_unfold,ens_id_rpkm,by="ENS_id") -> all_sample_sep_group_psm_unfold_add_rpkm
all_sample_sep_group_psm_unfold_add_rpkm %>% filter(Mean_rpkm>0) %>%  group_by(Group_id) %>%
  filter(n() == 1) %>% ungroup()-> psm_retained_group_based_on_transcript
# 基于此救回了50个微蛋白
psm_retained_group_based_on_transcript %>% count(Protein) %>% nrow()

# 4.合并结果
fwrite_c(psm_retained_group_based_on_transcript,"./psm_retained_group_based_on_transcript.txt")
## 1119
data_frame(
  ORF_id_trans=unique(c(psm_retained_group_based_on_transcript$Protein,in_house_res$ORF_id_trans))
  ) -> new_sep_list
nrow(new_sep_list)
fwrite_c(new_sep_list,"./new_sep_list.txt")
