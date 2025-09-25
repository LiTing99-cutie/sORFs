source("/home/user/data2/lit/bin/lit_utils.R")
source("/home/user/data3/lit/project/sORFs/sORFs.utils.R")
lib_text()
lib_plot()
setwd("/home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401")
# 读入进化分析的结果；注意这里只挑选了之前1070 list中类型为"ncORF","uORF","uoORF","dORF","doORF"的SEP拿去做进化分析
fread_c("/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/results_120/ancestors") -> ancestor
fread_c("/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/results_120/spec.out") -> spec.out
fread_c("/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/output/output/blastp_check/times_in_genome_trans.txt") -> times_in_genome_trans
colnames(times_in_genome_trans) <- c("orf_id","loc_times_blat","gene_times_blastn","gene_times_blastn_pc")
fread_c("/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/in_house_human_brain_denovo_check_20250410/output/output/blastp_check/outgroup_homolog.peptide_similarity.txt") -> outgroup_homolog.peptide_similarity
# table(ancestor$ev_age)
# table(ancestor$syn_age)
# table(spec.out$lineage)
# 读入newly detected SEP list
fread("./output/sep_add_basic_ms_ribo_info_group_retained.txt") -> all_sep
ancestor$orf_id <- sub(ancestor$orf_id,pattern = "__",replacement = ":")
spec.out$orf_id <- sub(spec.out$orf_id,pattern = "__",replacement = ":")
merge(ancestor[,c("orf_id","ev_age","syn_age","denovo")],spec.out[,c("orf_id","lineage")],all.y=T) -> m_evo
m_evo$denovo[is.na(m_evo$denovo)] <- 1
m_evo$ev_age[is.na(m_evo$ev_age)] <- "human"
m_evo$syn_age[is.na(m_evo$syn_age)] <- "human"
merge(m_evo,all_sep,by.x="orf_id",by.y="ORF_id_trans") -> all_sep_m_evo
fwrite_c(all_sep_m_evo,"./output/S8/all_sep_m_evo.txt")

merge(all_sep_m_evo,times_in_genome_trans,by="orf_id") %>% 
  merge(outgroup_homolog.peptide_similarity,by.x = "orf_id",by.y="Gene ID") -> all_sep_m_evo_add_homolog

filter(all_sep_m_evo_add_homolog,denovo==1) -> all_sep_m_evo_add_homolog_denovo

all_sep_m_evo_add_homolog_denovo %>% filter(loc_times_blat<=1 & gene_times_blastn<=1 & gene_times_blastn_pc<=1) -> res_filter_paralog

LEVELS <- c("human","hominoid","catarrhini","simiiformes","primates","primatomorpha","boreoeutheria",
            "placentalia","mammal")
res_filter_paralog$lineage <- factor(res_filter_paralog$lineage,levels=LEVELS)
res_filter_paralog$`Outgroup Homolog` <- factor(res_filter_paralog$`Outgroup Homolog`,levels=LEVELS)
res_filter_paralog$lineage_n <- as.numeric(res_filter_paralog$lineage)
res_filter_paralog$`Outgroup Homolog n` <- as.numeric(res_filter_paralog$`Outgroup Homolog`)
res_filter_paralog %>% filter(lineage_n>=`Outgroup Homolog n`) -> res_filter_paralog_ortholog
##### plot #####
### old result
all_sep_m_evo_add_homolog_denovo -> df
get_p <- function(df){
  "lineage" -> col
  data <- data.frame(table(df[,col]))
  colnames(data) <- c(col,"Freq")
  data$lineage <- factor(data$lineage,levels=LEVELS)
  p <- ggplot(data,aes_string(x = col, y = "Freq",fill=col)) + 
    geom_bar(stat = "identity",width = 0.8) +
    scale_y_continuous(expand = c(0, 0),limits=c(0,100)) +
    theme_3(rotate=T)+
    guides(fill = "none")
  p <- p+
    geom_text(aes(label = Freq), 
              position = position_dodge(0.8), # 确保标签对齐柱子
              vjust = -0.3,                  # 根据需要调整位置
              size = 5, color = "black")+
    coord_cartesian(clip = "off")
  p+scale_fill_discrete_divergingx() -> p
  return(p)
}
get_p(df) -> p
create_path("output/S8")
ggsave(p,filename = "./output/S8/denovo_lineage_specific.pdf",width = 6,height = 6)
### updated result
res_filter_paralog_ortholog -> df 
get_p(df)+scale_y_continuous(expand = c(0, 0),limits=c(0,30)) -> p
p
ggsave(p,filename = "./output/S8/denovo_lineage_specific_filter_homolog.pdf",width = 6,height = 6)
##### basic number #####
nrow(all_sep_m_evo_add_homolog)
nrow(all_sep_m_evo_add_homolog_denovo)
nrow(res_filter_paralog)
nrow(res_filter_paralog_ortholog)
res_filter_paralog_ortholog %>% filter(lineage %in% c("human","hominoid","catarrhini","simiiformes","primates","primatomorpha")) %>% nrow()
res_filter_paralog_ortholog %>% filter(lineage %in% c("human","hominoid")) %>% nrow()
##### hominoid-specific genes #####
res_filter_paralog_ortholog %>% filter(lineage %in% c("human","hominoid")) -> denovo_genes
table(denovo_genes$ORF_type_3) %>% sort()
table(denovo_genes$Scodon) %>% sort()
