# 将riborf的ID名整理成统一格式；去掉不同转录本编码的相同ORF，选出代表性的转录本
# 1.加载R包
args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R")
lib_text()
lib_plot()

# 2.指定路径
# 折叠后的结果
collase_trans_orf_gpe_path <- "./annot/RiboORF/mm39/candidateORF.genepred.callapse.stop_codon.txt"
# 所有的orf的蛋白质序列
trans_orf_fa_path <- "./annot/RiboORF/mm39/candidateORF.prot.tab"
output_path <- "./annot/RiboORF/mm39/"

trans_mul_feature_path <- "/home/user/data3/lit/resource/gtf/mouse/mm39/Ensembl_106_Gencode_vM29_Mouse_Transcript_stable_ID_version_TSL_APPRIS_Transcript_length_Gene_stable_ID_version_Gene_name_CDS_length_Transcript_type.txt"
script_path <- "./S3.1b1.Run.Organize_transcript_meta.R"
uniprot_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/uniprot/mouse/uniprotkb_Mus_musculus_reviewed_canonical_and_isoform.rmdup.seq"
ncbi_seq_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/annot/NCBI_refseq/mm39/GCF_000001635.27_GRCm39_protein.rmdup.seq"

# 3.读取文件并合并
# # 数据框比较大，因此先读取1000行做测试
# fread(collase_trans_orf_gpe_path, nrows = 1000,data.table = F) -> collase_trans_orf_gpe
# fread(trans_sorf_fa_path, nrows = 1000,data.table = F,sep = '\t') %>% .[,c(1,2)]-> trans_sorf_fa

# 读取两个文件并合并
## 读取数据时尽量只选择必要的列
fread(collase_trans_orf_gpe_path,select = c("V1", "V6", "V7"),data.table = F) -> collase_trans_orf_gpe
collase_trans_orf_gpe %>% rename("RibORF_id"="V1","codon5"="V6","codon3"="V7") -> collase_trans_orf_gpe
fread(trans_orf_fa_path,data.table = F,sep = '\t') %>% .[,c(1,2)] -> trans_orf_fa
colnames(trans_orf_fa) <- c("RibORF_id","Seq") 
merge(trans_orf_fa,collase_trans_orf_gpe,by="RibORF_id") -> collapse_trans_orf

# 4.Uniform ID名称
# 数据框比较大，因此先读取1000行做测试
# head(collapse_trans_orf,100) -> df

# 替换ID名称
tf_riborf_id <- function(df){
  # tf为transform的意思，基于riborf的结果生成uniform的两个id
  # 输出的数据框必须有RibORF_id以及seq序列，codon5以及codon3的位置
  df <- df %>%
    separate(col = "RibORF_id", into = c("ENS_id", "Chr", "Strand", "Rank", "Span","RelaStart", "RelaEnd", "Type", "Scodon"), sep = "[|:]", remove = FALSE)
  ORF_id_trans <- paste0(df$ENS_id,df$Strand,df$Chr,
                   ":",df$codon5,"-",df$codon3)
  df$ORF_id_trans <- ORF_id_trans
  df$Location <- str_extract(df$ORF_id_trans,"[+-]chr\\w+:\\d+[+-]\\d+")
  df$ORF_id_seq <- paste0(df$Location,":",df$Seq)
  return(df)
}

tf_riborf_id(collapse_trans_orf) -> all_orfs

# 5.过滤结果
top_trans_filter <- function(all_orfs){
  # 选择代表性转录本，并且去除和uniprot完全overlap的蛋白质
  # 输入all_orfs，需要有ORF_id_seq以及ENS_id列
  
  # ORF_id可能对应多个ENS_id
  distinct(all_orfs,ORF_id_seq,ENS_id) -> ORF_id_ENS_id
  # 3.整合从biomart中下载的转录本的信息并与sORF合并
  source(script_path)
  trans_mul_feature_add <- fread(trans_mul_feature_path)
  merge(ORF_id_ENS_id,trans_mul_feature_add[,c("Transcript_ID","TSL","APPRIS","CDS_l")],by.x = "ENS_id",by.y = "Transcript_ID") -> ORF_id_ENS_id_merge_feature
  
  # 4.筛选出证据最强的转录本
  # 设置TSL和APPRIS的等级排序
  tsl_levels <- c("tsl1", "tsl2", "tsl3", "tsl4", "tsl5", "tslNA")
  appris_levels <- c("principal1", "principal2", "principal3","principal4", "principal5", "alternative1", "alternative2", "")
  # 转换TSL和APPRIS为因子并指定等级
  ORF_id_ENS_id_merge_feature <- ORF_id_ENS_id_merge_feature %>%
    mutate(TSL = factor(TSL, levels = tsl_levels),
           APPRIS = factor(APPRIS, levels = appris_levels))
  # 进行排序并筛选出证据最强的转录本
  ORF_id_ENS_id_top_trans <- ORF_id_ENS_id_merge_feature %>%
    group_by(ORF_id_seq) %>%
    arrange(APPRIS,TSL,desc(CDS_l)) %>%  # 按TSL、APPRIS等级以及cds_length降序排序
    slice(1) %>%  # 保留每组的第一行
    ungroup()     # 解除分组
  
  # 5.得到基于转录本去重的所有sORFs
  merge(all_orfs,ORF_id_ENS_id_top_trans[,c("ENS_id","ORF_id_seq")],c("ENS_id","ORF_id_seq")) -> all_orfs_top_trans
  
  # 6.去掉和NCBI或者uniprot蛋白质序列相同的sORF
  ncbi_seq <- fread(ncbi_seq_path,header = F)
  uniprot_seq <- fread(uniprot_seq_path,header = F)
  all_orfs_top_trans$Seq %in% c(ncbi_seq$V1,uniprot_seq$V1) -> idx
  all_orfs_top_trans[!idx,] -> all_orfs_filtered
  
  # 返回经过过滤的所有orf
  return(all_orfs_filtered)
  
}

top_trans_filter(all_orfs) -> all_orfs_top_trans_filter

all_orfs_top_trans_filter$Length <- str_length(all_orfs_top_trans_filter$Seq)

all_orfs_top_trans_filter %>% 
  .[,c("ORF_id_trans","ORF_id_seq","Seq","Length","Chr","codon5","codon3","Strand","ENS_id","Type","Scodon")] -> all_orfs_top_trans_filter

# 6.过滤得到小肽
summary(all_orfs_top_trans_filter)
all_orfs_top_trans_filter %>% filter(Length>=6 & Length<=150) -> trans_based_sorfs
nrow(all_orfs_top_trans_filter)-nrow(trans_based_sorfs)
fwrite_c(trans_based_sorfs,paste0(output_path,"/","trans_based_sorfs.txt"))
