source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
# 将riborf的ID名整理成统一格式
setwd("/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/trans_based_database")
trans_orf_fa_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/mouse_run_add_assemble_20250125/annotation/RibORF_annot/candidateORF.prot.tab"
trans_orf_gpe_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/mouse_run_add_assemble_20250125/annotation/RibORF_annot/candidateORF.genepred.txt"
# 1. 加载所需数据
trans_orf_gpe <- fread(trans_orf_gpe_path, select = c("V1", "V6", "V7"))
setnames(trans_orf_gpe, c("V1", "V6", "V7"), c("RibORF_id", "codon5", "codon3"))
trans_orf_fa <- fread(trans_orf_fa_path, select = c(1, 2),sep = '\t')
setnames(trans_orf_fa, c("RibORF_id", "Seq"))
all_orfs <- merge(trans_orf_fa, trans_orf_gpe, by = "RibORF_id")
# 2. 统一ID名称
tf_riborf_id <- function(df) {
  # 拆分RibORF_id生成统一ID
  df[, c("ENS_id", "Chr", "Strand", "Rank", "Span", "RelaStart", "RelaEnd", "Type", "Scodon") :=
       tstrsplit(RibORF_id, split = "[|:]", fixed = FALSE)]
  df[, ORF_id_trans := paste0(ENS_id, Strand, Chr, ":", codon5, "-", codon3)]
  df[, Location := str_extract(ORF_id_trans,"[+-]chr\\w+:\\d+[+-]\\d+")]
  df[, ORF_id_seq := paste0(Location, ":", Seq)]
  return(df)
}
all_orfs_tf <- tf_riborf_id(all_orfs)
saveRDS(all_orfs_tf,"./output/all_orfs_tf.rds")
# 3.Check是否所有的ribo orfs在基于trans的database中
ribo_sorfs_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/mouse_brain_output_20241011/Ribo_ORFs_add_assemble_20250125/filter/nonCano.sorf.filtered.add_meta.orf_type_annotated.txt"
fread_c(ribo_sorfs_path) -> ribo_sorfs
table(ribo_sorfs$ORF_id_trans %in% all_orfs_tf$ORF_id_trans)
table(ribo_sorfs$ORF_id_seq %in% all_orfs_tf$ORF_id_seq)


