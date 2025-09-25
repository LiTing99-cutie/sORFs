
# 1.加载R包
# args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()

trans_based_sorfs_path <- "annot/sORF_database/Transcriptome_based/trans_based_sorfs.txt"
ribo_base_sorfs_path <- "mouse_brain_output_20241011/Ribo-ORFs-20241218/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt"

trans_based_sorfs <- fread(trans_based_sorfs_path)
ribo_based_sorfs <- fread(ribo_base_sorfs_path)

# 2.比较基于Ribo-seq和基于转录组得到的两个库的区别
ribo_based_sorfs$ORF_id_seq %in% trans_based_sorfs$ORF_id_seq %>% sum()

venn_plot_n(list(
  trans_based_sorfs=trans_based_sorfs$ORF_id_seq,
  ribo_based_sorfs=ribo_based_sorfs$ORF_id_seq
))

# 3.得到只在Trans-based库中得到的小肽
trans_based_sorfs[!trans_based_sorfs$ORF_id_seq %in% ribo_based_sorfs$ORF_id_seq] -> trans_based_specific_sorfs
fwrite_c(trans_based_specific_sorfs,paste0(output_path,"/","trans_based_specific_sorfs.txt"))
