library(profvis)
profvis({
  collase_trans_orf_gpe <- fread(collase_trans_orf_gpe_path, select = c("V1", "V6", "V7"))
  setnames(collase_trans_orf_gpe, c("V1", "V6", "V7"), c("RibORF_id", "codon5", "codon3"))
  
  trans_orf_fa <- fread(trans_orf_fa_path, select = c(1, 2),sep = '\t')
  setnames(trans_orf_fa, c("RibORF_id", "Seq"))
  
  # collapse_trans_orf <- data.table::merge(trans_orf_fa, collase_trans_orf_gpe, by = "RibORF_id")
})

p <- profvis({
  collapse_trans_orf <- data.table::merge.data.table(trans_orf_fa, collase_trans_orf_gpe, by = "RibORF_id")
})

p_1 <- profvis({
  fread(collase_trans_orf_gpe_path,select = c("V1", "V6", "V7"),data.table = F) -> collase_trans_orf_gpe
  collase_trans_orf_gpe %>% rename("RibORF_id"="V1","codon5"="V6","codon3"="V7") -> collase_trans_orf_gpe
  fread(trans_orf_fa_path,data.table = F,sep = '\t') %>% .[,c(1,2)] -> trans_orf_fa
  colnames(trans_orf_fa) <- c("RibORF_id","Seq") 
})

system.time(
  {
    fread(collase_trans_orf_gpe_path,select = c("V1", "V6", "V7"),data.table = F) -> collase_trans_orf_gpe
    collase_trans_orf_gpe %>% rename("RibORF_id"="V1","codon5"="V6","codon3"="V7") -> collase_trans_orf_gpe
    fread(trans_orf_fa_path,data.table = F,sep = '\t') %>% .[,c(1,2)] -> trans_orf_fa
    colnames(trans_orf_fa) <- c("RibORF_id","Seq") 
  }
)

system.time(
  {
    collase_trans_orf_gpe <- fread(collase_trans_orf_gpe_path, select = c("V1", "V6", "V7"))
    setnames(collase_trans_orf_gpe, c("V1", "V6", "V7"), c("RibORF_id", "codon5", "codon3"))
    
    trans_orf_fa <- fread(trans_orf_fa_path, select = c(1, 2),sep = '\t')
    setnames(trans_orf_fa, c("RibORF_id", "Seq"))
  }
)

system.time(
{tf_riborf_id <- function(df) {
  # 拆分RibORF_id生成统一ID
  df[, c("ENS_id", "Chr", "Strand", "Rank", "Span", "RelaStart", "RelaEnd", "Type", "Scodon") :=
       tstrsplit(RibORF_id, split = "[|:]", fixed = FALSE)]
  df[, ORF_id_trans := paste0(ENS_id, Strand, Chr, ":", codon5, "-", codon3)]
  df[, Location := str_extract(ORF_id_trans,"[+-]chr\\w+:\\d+[+-]\\d+")]
  df[, ORF_id_seq := paste0(Location, ":", Seq)]
  return(df)
}

all_orfs <- tf_riborf_id(collapse_trans_orf)}
)

system.time(
  {tf_riborf_id <- function(df){
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
  
  tf_riborf_id(collapse_trans_orf) -> all_orfs}
)

system.time({
# 3. 过滤结果
top_trans_filter <- function(all_orfs) {
  # ORF_id可能对应多个ENS_id
  ORF_id_ENS_id <- unique(all_orfs[, .(ORF_id_seq, ENS_id)])
  
  # 读取转录本特征信息并合并
  trans_mul_feature_add <- fread(trans_mul_feature_path)
  ORF_id_ENS_id_merge_feature <- merge(
    ORF_id_ENS_id,
    trans_mul_feature_add[, .(Transcript_ID, TSL, APPRIS, CDS_l)],
    by.x = "ENS_id",
    by.y = "Transcript_ID"
  )
  
  # 设置TSL和APPRIS的等级排序
  tsl_levels <- c("tsl1", "tsl2", "tsl3", "tsl4", "tsl5", "tslNA")
  appris_levels <- c("principal1", "principal2", "principal3", "principal4", "principal5", "alternative1", "alternative2", "")
  
  ORF_id_ENS_id_merge_feature[, TSL := factor(TSL, levels = tsl_levels)]
  ORF_id_ENS_id_merge_feature[, APPRIS := factor(APPRIS, levels = appris_levels)]
  
  # 按条件排序并筛选出证据最强的转录本
  ORF_id_ENS_id_top_trans <- ORF_id_ENS_id_merge_feature[
    order(ORF_id_seq, APPRIS, TSL, -CDS_l),
    .SD[1],
    by = ORF_id_seq
  ]
  
  # 得到基于转录本去重的所有sORFs
  all_orfs_top_trans <- merge(
    all_orfs,
    ORF_id_ENS_id_top_trans[, .(ENS_id, ORF_id_seq)],
    by = c("ENS_id", "ORF_id_seq")
  )
  
  # 去掉和NCBI或Uniprot蛋白质序列相同的sORF
  ncbi_seq <- fread(ncbi_seq_path, header = FALSE)
  uniprot_seq <- fread(uniprot_seq_path, header = FALSE)
  
  all_orfs_top_trans <- all_orfs_top_trans[!Seq %in% c(ncbi_seq$V1, uniprot_seq$V1)]
  
  return(all_orfs_top_trans)
}

# 过滤并保存结果
all_orfs_top_trans_filter <- top_trans_filter(all_orfs)
})
