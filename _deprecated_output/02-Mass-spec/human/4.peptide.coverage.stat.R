source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
setwd("/home/user/data3/lit/project/sORFs/02-Mass-spec/human/")
novel_sep_add_info_path <- "/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/sep_add_basic_ms_ribo_info_group_retained.txt"
##### 肽段在SEP上相对位置的分布【统计】 ##### 
fread_c("./psm_retained_group_based_on_transcript.txt") -> psm
fread_c("./sep_unique_psm.txt") -> psm_1
cols <- c("Protein","Peptide")
rbind(psm[,cols],psm_1[,cols]) %>% distinct(Protein,Peptide) -> pro_pep
fread(novel_sep_add_info_path) -> novel_sep
# 需要更新ID
pro_pep %>% dplyr::rename(ORF_id_trans=Protein) -> tmp
update_orf_id_trans(tmp) -> pro_pep_1
# 获取起始密码子等信息
merge(pro_pep_1,novel_sep,by="ORF_id_trans") -> all_sep_ms_pep
# 创建输出路径
create_path("./S4")
# 写出ID转换后的novel SEP的肽段支持文件
fwrite_c(all_sep_ms_pep,"./S4/novel_sep_ms_pep.txt")