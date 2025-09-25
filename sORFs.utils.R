#!/bin/env Rscript

################################################
#File Name: sORFs.utils.R
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Tue 01 Apr 2025 07:54:38 PM CST
################################################


##### 读取数据 ######
# 从蛋白质的结果中去掉去掉污染蛋白
fil_contam <- function(df=df,species="mouse"){
  # 读取污染蛋白ID
  fread("/rd1/user/lit/project/sORFs/custom_database/crap.Entry.Name.txt") -> contam_id
  df$`Entry Name` %in% contam_id -> idx
  df[!idx,] -> df_fil_contam
  # 有些蛋白也是污染蛋白，但是并没有完全包括在contam_id中，得到物种为小鼠或者""也就是新鉴定小肽ID对应的蛋白
  if(species=="mouse"){
    df_fil_contam %>% filter(.,Organism=="Mus musculus" | Organism=="" | is.na(Organism)) -> df_fil_contam
  }
  return(df_fil_contam)
}

read_file_fil_contam <- function(path=path,species="mouse"){
  # 读取蛋白质结果，并从蛋白质的结果中去掉去掉污染蛋白，并从路径中提取样本名
  df <- fread(path,data.table = F)
  if(nrow(df)>1){
    if(species=="mouse"){
      df <- fil_contam(df,"mouse")
    }else{
      df <- fil_contam(df)
    }
    df$Sample <- str_extract(path, "(?<=min_)\\w+(?=_Slot)")
  }
  return(df)
}

read_file <- function(path=path){
  # 读取结果，并从路径中提取样本名
  df <- fread(path,data.table = F)
  if(nrow(df)>1){
    df$Sample <- str_extract(path, "(?<=min_)\\w+(?=_Slot)")
  }
  return(df)
}

##### 筛选 #####
# 得到小肽对应的数据框
get_sep <- function(df){
  filter(df,Length<=150)
}

# 得到新鉴定小肽对应的数据框
get_new_sep <- function(df){
  return(df[grepl("^ENS|^STRG",df$Protein),])
}

get_uni <- function(df){
  filter(df,`Unique Peptides`>0) -> df
  return(df)
}

##### 合并和基本统计 ######
# 从路径中combine所有的蛋白质结果，并从蛋白质的结果中去掉去掉污染蛋白
get_combined_protein_df_f_path <- function(path){
  protein_files <- list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T)
  lapply(protein_files,read_file_fil_contam,species="mouse") -> l
  names(l) <- str_extract(protein_files, "(?<=min_)\\w+(?=_Slot)")
  merged_protein <- do.call(rbind, l) 
  return(merged_protein)
}

read_psm <- function(psm_path) {
  # 检查文件是否只有标题行
  if (length(readLines(psm_path)) == 1) {
    return(NULL)  # 返回空数据框以便后续合并
  }else{
  df <- read.table(psm_path, header = TRUE, sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE)
  df$Sample <- str_extract(psm_path, "(?<=min_)\\w+(?=_Slot)")
  return(df)
  }
}

get_total_psm <- function(path) {
  # 得到某个路径下所有的psm结果，并且append起来
  files <- list.files(path, pattern = "^psm.tsv", full.names = TRUE,recursive = T)
  files %>% lapply(.,read_psm) %>%
    # 过滤掉list中的null值
    Filter(Negate(is.null), .) %>% 
    do.call(rbind,.) -> total_psm
  total_psm
}


##### 可视化 ######
plot_certain_sample_venn <- function(merged_df=merged_df,compare_lst=compare_lst){
  sample_split_protein <- split(merged_df[, "Protein", drop = FALSE], merged_df$Sample)
  lapply(sample_split_protein, function(x){return(x$Protein)}) -> sample_split_protein_chr
  venn_plot_n(sample_split_protein_chr[compare_lst]) -> p1
  return(p1)
}

##### 小肽统计相关 ######
get_grouped_protein <- function(df){
  # 展开Indistinguishable Proteins中的蛋白
  # 输入：protein.tsv数据框，其中Indistinguishable Proteins这一列非空
  # 输出：展开了indistinguishable Proteins的数据框，其中蛋白质名字被替换，且subset new_sep出来，但是Length长度需要进一步修改
  df %>% separate_rows(`Indistinguishable Proteins`, convert = TRUE,sep=", ") -> tmp
  tmp %>% mutate(Protein=`Indistinguishable Proteins`,
                 `Protein ID`=`Indistinguishable Proteins`,
                 `Entry Name`=`Indistinguishable Proteins`,
                 `Protein Description`=`Indistinguishable Proteins`) -> tmp_1
  return(tmp_1[grepl("^ENS|^STRG",tmp_1$Protein),])
}

get_leveled_new_sep_df <- function(protein_file){
  # 对于每个蛋白质结果文件，得到一个new_sep的df
  # 输入：protein.tsv
  # 输出：样本中新鉴定到的所有小肽，展开了Indistinguishable Proteins，并且Confidence_Level中标注了可信度
  read_file(protein_file) -> protein
  # 得到新鉴定的小肽
  get_new_sep(protein) -> s_pep
  # 对新鉴定的小肽进行分组，分为被unique peptide支持，protein group中的leader protein，protein group中只有new_sep，protein group中两者兼有或者只有注释蛋白
  mutate(
    s_pep,
    sp_present = grepl("sp", s_pep$`Indistinguishable Proteins`),
    starts_with_ENS = grepl("^ENS|^STRG", s_pep$`Indistinguishable Proteins`),
    contains_ENS = grepl(" ENS| STRG", s_pep$`Indistinguishable Proteins`),
    Group_type = case_when(
      `Unique Peptides`>0 ~ "Unique",
      `Unique Peptides`==0 & `Indistinguishable Proteins`=="" ~ "Group_none",
      !sp_present & (starts_with_ENS | contains_ENS) ~ "Group_only_new_sep",
      sp_present ~ "Group_mixed"
    )
  )  %>% 
    select(-sp_present, -starts_with_ENS, -contains_ENS) -> s_pep_add_type
  # 得到不同证据支持程度的小肽
  s_pep_add_type %>% filter(Group_type=="Unique") %>% mutate(Confidence_Level=1) -> Level_1_new_sep
  s_pep_add_type %>% filter(Group_type=="Group_none") %>% mutate(Confidence_Level=2) -> Level_2_new_sep
  ## 对于protein列和Indistinguishable Proteins同时含有新鉴定小肽的情况，需要把原数据和展开后的数据进行合并
  s_pep_add_type %>% filter(Group_type=="Group_only_new_sep") -> df 
  rbind(df,get_grouped_protein(df)) %>% mutate(Confidence_Level=3) -> Level_3_new_sep
  s_pep_add_type %>% filter(Group_type=="Group_mixed") -> df_1
  ## 对于protein列为注释蛋白，Indistinguishable Proteins中可能含有新鉴定小肽，展开
  protein %>% filter(!grepl("^ENS|^STRG", protein$Protein) & `Indistinguishable Proteins`!="") -> df_2
  rbind(df_1,get_grouped_protein(df_1),get_grouped_protein(df_2)) %>% 
    mutate(Confidence_Level=4) -> Level_4_new_sep
  rbind(Level_1_new_sep,Level_2_new_sep,Level_3_new_sep,Level_4_new_sep) -> All_level_sep  
  All_level_sep %>% select(-Group_type) -> All_level_sep
  return(All_level_sep)
}

get_sample_ms_stat_from_path <- function(PATH){
  # 输入：结果路径，子路径中包含了protein.tsv、peptide.tsv、ion.tsv以及psm.tsv，可以输入多个路径构成的向量
  # 输出：每个样本中的质谱统计数据：蛋白总数，肽段总数，以及肽段谱图匹配PSM总数
  protein_files <- c()
  peptide_files <- c()
  ion_files <- c()
  psm_files <- c()
  for (path in PATH){
    protein_files <- c(list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T),protein_files)
    peptide_files <- c(list.files(path, pattern = "^peptide.tsv", full.names = TRUE,recursive = T),peptide_files)
    psm_files <- c(list.files(path, pattern = "^psm.tsv", full.names = TRUE,recursive = T),psm_files)
  }
  df <- data.frame(
    Sample=str_extract(protein_files, "(?<=min_)\\w+(?=_Slot)"),
    Protein_Group_N=sapply(protein_files, function(x){read_file_fil_contam(x) %>% nrow()}),
    Peptdide_N=sapply(peptide_files, function(x){read_file(x) %>% nrow()}),
    PSM_N=sapply(psm_files, function(x){read.table(x, header = TRUE, sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE) %>% nrow()})
  )
  # 原始的行名为文件名，去除
  rownames(df) <- NULL
  return(df)
}

get_stat_from_sample <- function(PATH){
  # 输入：结果路径，子路径中包含了protein.tsv，可以输入多个路径构成的向量
  # 输出：小肽总体统计；分样本小肽以及质谱数量统计；小肽构成的数据框
  protein_files <- c()  
  for (path in PATH){
    protein_files <- c(list.files(path, pattern = "^protein.tsv", full.names = TRUE,recursive = T),protein_files)}
  lapply(protein_files, get_leveled_new_sep_df) %>% do.call(rbind,.) -> all_sample_all_level_sep
  # 统计小肽结果
  ## 统计总数
  # 统计不同置信度的小肽时，如果在某个样本中置信度是1，但是在其他样本中置信度是2，那往高了选取，在任意一个样本中置信度是1，那么就是1
  all_sample_all_level_sep %>% group_by(Protein) %>% filter(Confidence_Level==min(Confidence_Level)) %>% 
    distinct(Protein,.keep_all = T) %>% ungroup() %>% count(Confidence_Level) -> tmp
  all_sep_sta <- data.frame(Feature=c("New_sep_level_1_N","New_sep_level_2_N",
                                      "New_sep_level_3_N","New_sep_level_4_N","New_sep_total_N"),
                            N=c(tmp$n,n_distinct(all_sample_all_level_sep$Protein)))
  all_sep_sta
  ## 分样本进行统计
  all_sample_all_level_sep %>% count(Sample,Confidence_Level) %>% dcast(Sample~Confidence_Level) -> sample_sep_sta
  sample_sep_sta[is.na(sample_sep_sta)] <- 0
  mutate(sample_sep_sta,New_sep_total_N=rowSums(sample_sep_sta[2:5])) -> sample_sep_sta
  colnames(sample_sep_sta) <- c("Sample","New_sep_level_1_N","New_sep_level_2_N",
                                "New_sep_level_3_N","New_sep_level_4_N","New_sep_total_N")
  sample_sep_sta
  # 统计质谱搜库结果
  ## 分样本进行统计
  get_sample_ms_stat_from_path(PATH) -> sample_ms_sta
  merge(sample_sep_sta,sample_ms_sta,by="Sample") -> sample_sep_ms_sta  
  
  return(list(sample_sep_ms_sta=sample_sep_ms_sta,all_sep_sta=all_sep_sta,all_sample_all_level_sep=all_sample_all_level_sep))
}

# 绘制累积曲线
cumu_sep_plot <- function(df){
  # 统计累积 unique 蛋白质数量
  # 输入一个蛋白质和样本对应的数据框
  df$Sample <- factor(df$Sample,levels = sample_order)
  arrange(df,Sample) -> df
  df_1 <- df %>% distinct(Protein,.keep_all = T)
  cumulative_data <- df_1 %>%
    group_by(Sample) %>%
    summarise(Proteins = list(unique(Protein))) %>%
    mutate(CumulativeUniqueProteins = cumsum(sapply(Proteins, length)),
           SampleCount = row_number())
  
  # 绘制累积频数曲线
  ggplot(cumulative_data, aes(x = SampleCount, y = CumulativeUniqueProteins)) +
    geom_line() +
    geom_point() +
    theme_3() -> p_1
  
  # 确保 Sample 顺序
  cumulative_data$Sample <- factor(cumulative_data$Sample, levels = unique(cumulative_data$Sample))
  # cumulative_data$Sample <- factor(cumulative_data$Sample, levels = sample_order)
  ggplot(cumulative_data, aes(x = Sample, y = CumulativeUniqueProteins,group=1)) +
    geom_line() +
    geom_point() +
    theme_3(rotate = T) -> p_2
  
  return(list(cumulative_data=cumulative_data,p_1=p_1,p_2=p_2))
}
cumu_sep_plot_v1 <- function(df,sample_order){
  # 统计累积 unique 蛋白质数量
  # 输入一个蛋白质和样本对应的数据框
  df$Sample <- factor(df$Sample,levels = sample_order)
  arrange(df,Sample) -> df
  df_1 <- df %>% distinct(Protein,.keep_all = T)
  cumulative_data <- df_1 %>%
    group_by(Sample) %>%
    summarise(Proteins = list(unique(Protein))) %>%
    mutate(CumulativeUniqueProteins = cumsum(sapply(Proteins, length)),
           SampleCount = row_number())
  
  # 绘制累积频数曲线
  ggplot(cumulative_data, aes(x = SampleCount, y = CumulativeUniqueProteins)) +
    geom_line() +
    geom_point() +
    theme_3() -> p_1
  
  # 确保 Sample 顺序
  cumulative_data$Sample <- factor(cumulative_data$Sample, levels = unique(cumulative_data$Sample))
  ggplot(cumulative_data, aes(x = Sample, y = CumulativeUniqueProteins,group=1)) +
    geom_line() +
    geom_point() +
    theme_3(rotate = T) -> p_2
  
  return(list(cumulative_data=cumulative_data,p_1=p_1,p_2=p_2))
}
cumu_sep_peptide_plot <- function(df,sample_order){
  # 统计累积 unique 肽段数量
  # 输入一个肽段和样本对应的数据框
  df$Sample <- factor(df$Sample,levels = sample_order)
  arrange(df,Sample) -> df
  df_1 <- df %>% distinct(Peptide,.keep_all = T)
  cumulative_data <- df_1 %>%
    group_by(Sample) %>%
    summarise(Peptides = list(unique(Peptide))) %>%
    mutate(CumulativeUniquePeptides = cumsum(sapply(Peptides, length)),
           SampleCount = row_number())
  
  # 绘制累积频数曲线
  ggplot(cumulative_data, aes(x = SampleCount, y = CumulativeUniquePeptides)) +
    geom_line() +
    geom_point() +
    theme_3() -> p_1
  
  # 确保 Sample 顺序
  cumulative_data$Sample <- factor(cumulative_data$Sample, levels = unique(cumulative_data$Sample))
  ggplot(cumulative_data, aes(x = Sample, y = CumulativeUniquePeptides,group=1)) +
    geom_line() +
    geom_point() +
    theme_3(rotate = T) -> p_2
  
  return(list(cumulative_data=cumulative_data,p_1=p_1,p_2=p_2))
}
cumu_combo_plot <- function(df, sample_order, target_col = "Protein", 
                            show_bar = TRUE, 
                            bar_color = hcl.colors(3, "Dark 3")[3],
                            line_color = hcl.colors(3, "Dark 3")[1]) {
  
  # 检查必要列
  if (!all(c("Sample", target_col) %in% colnames(df))) {
    stop(paste("输入数据必须包含Sample和", target_col, "列"))
  }
  
  # 设置样本顺序并去重
  df$Sample <- factor(df$Sample, levels = sample_order)
  df <- arrange(df, Sample)
  df_unique <- df %>% distinct(!!sym(target_col), .keep_all = TRUE)
  
  # 计算累积数据
  ## 捞一下没有unique蛋白质的样本
  tmp_all_sample <- data.frame(unique(df$Sample))
  colnames(tmp_all_sample) <-  "Sample"
  df_unique %>%
  group_by(Sample) %>%
  summarise(Targets = list(unique(!!sym(target_col)))) -> tmp
  merge(tmp,tmp_all_sample,by="Sample",all = T) -> tmp_1
  tmp_1$Sample <- factor(tmp_1$Sample, levels = sample_order)
  tmp_1 <- arrange(tmp_1, Sample)
  tmp_1 %>%  mutate(
      CumulativeUnique = cumsum(sapply(Targets, length)),
      SampleCount = row_number()
    ) -> cumulative_data
  
  # 计算每样本独立计数（用于柱状图）
  sample_counts <- df %>% 
    distinct(Sample, !!sym(target_col)) %>% 
    count(Sample)

  # 合并数据
  plot_data <- left_join(cumulative_data, sample_counts, by = "Sample")
  
  # 基础绘图
  p <- ggplot(plot_data, aes(x = Sample)) +
    geom_line(aes(y = CumulativeUnique, group = 1), 
              size = 1, color = line_color) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_3(rotate = TRUE) +
    labs(y = paste("Cumulative Unique", target_col))
  
  # 可选柱状图
  if (show_bar) {
    p <- p + geom_bar(aes(y = n), stat = "identity", 
                      fill = bar_color, width = 0.8)
  }
  
  # 返回结果
  return(list(
    plot_data = plot_data,
    version_sampleName = p,
    version_sampleNum = p %+% aes(x = SampleCount)  # 切换为样本序号版本
  ))
}

##### 导出谱图 ######
reformat_spectrum <- function(psm){
  # 增加一列
  psm$Spectrum %>% stringr::str_split_fixed(.,"\\.",4) %>% .[,2] %>% gsub("^0+", "", .) -> scan_number
  psm$Spectrum %>% stringr::str_split_fixed(.,"\\.",4) %>% .[,4] -> charge
  paste0(stringr::str_split_fixed(psm$Spectrum,"\\.",4) %>% .[,1],".",
         scan_number,".",scan_number,".",charge) -> spectrum
  psm$Spectrum_1 <- spectrum
  return(psm)
}

###### 小肽类型assign #####
get_orf_type_1_v0 <- function(Trans_type,Strand,Start,End,CDS_start,CDS_end){
  if(Trans_type!="protein_coding"){
    return("ncORF")}
  else {
    if (Strand == "+"){
      if (Start == CDS_start && End == CDS_end){
        return("Canonical")
      }else if (End <= CDS_start){
        return("uORF")
      }else if (Start < CDS_start && End > CDS_start && End < CDS_end){
        return("uoORF")
      }else if (Start < CDS_start && End == CDS_end){
        return("Extension")
      }else if (Start == CDS_start && End > CDS_end){
        return("Readthrough")
      }else if (Start < CDS_start && End > CDS_end){
        return("Containing")
      }else if (Start == CDS_start){
        return("Other_CDS_variant")
      }else if (Start >= CDS_end){
        return("dORF")
      }else if (End > CDS_end && Start > CDS_start && Start < CDS_end){
        return("doORF")
      }else if (Start > CDS_start && End < CDS_end){
        return("Internal")
      }else if (Start > CDS_start && End == CDS_end){
        return("Truncation")
      }
    }
    else if (Strand == "-"){
      if (Start == CDS_start && End == CDS_end){
        return("Canonical")
      }else if (Start >= CDS_end){
        return("uORF")
      }else if (End > CDS_end && Start > CDS_start && Start < CDS_end){
        return("uoORF")
      }else if (End > CDS_end && Start == CDS_start){
        return("Extension")
      }else if (End == CDS_end && Start < CDS_start){
        return("Readthrough")
      }else if (End > CDS_end && Start < CDS_start){
        return("Containing")
      }else if (End == CDS_end){
        return("Other_CDS_variant")
      }else if (End <= CDS_start){
        return("dORF")
      }else if (Start < CDS_start && End > CDS_start && End < CDS_end){
        return("doORF")
      }else if (Start > CDS_start && End < CDS_end){
        return("Internal")
      }else if (Start == CDS_start && End < CDS_end){
        return("Truncation")
      }
    }
  }
  return("other")
}
get_orf_type_1 <- function(Trans_type,Strand,Start,End,CDS_start,CDS_end){
  if(CDS_start==CDS_end){
    return("ncORF")}
  else {
    if (Strand == "+"){
      if (Start == CDS_start && End == CDS_end){
        return("Canonical")
      }else if (End <= CDS_start){
        return("uORF")
      }else if (Start < CDS_start && End > CDS_start && End < CDS_end){
        return("uoORF")
      }else if (Start < CDS_start && End == CDS_end){
        return("Extension")
      }else if (Start == CDS_start && End > CDS_end){
        return("Readthrough")
      }else if (Start < CDS_start && End > CDS_end){
        return("Containing")
      }else if (Start == CDS_start){
        return("Other_CDS_variant")
      }else if (Start >= CDS_end){
        return("dORF")
      }else if (End > CDS_end && Start > CDS_start && Start < CDS_end){
        return("doORF")
      }else if (Start > CDS_start && End < CDS_end){
        return("Internal")
      }else if (Start > CDS_start && End == CDS_end){
        return("Truncation")
      }
    }
    else if (Strand == "-"){
      if (Start == CDS_start && End == CDS_end){
        return("Canonical")
      }else if (Start >= CDS_end){
        return("uORF")
      }else if (End > CDS_end && Start > CDS_start && Start < CDS_end){
        return("uoORF")
      }else if (End > CDS_end && Start == CDS_start){
        return("Extension")
      }else if (End == CDS_end && Start < CDS_start){
        return("Readthrough")
      }else if (End > CDS_end && Start < CDS_start){
        return("Containing")
      }else if (End == CDS_end){
        return("Other_CDS_variant")
      }else if (End <= CDS_start){
        return("dORF")
      }else if (Start < CDS_start && End > CDS_start && End < CDS_end){
        return("doORF")
      }else if (Start > CDS_start && End < CDS_end){
        return("Internal")
      }else if (Start == CDS_start && End < CDS_end){
        return("Truncation")
      }
    }
  }
  return("other")
}
get_orf_type_2 <- function(sORF){
  # 生成ORF_type
  mapply(get_orf_type_1,sORF$Transcript_type,sORF$Strand,sORF$Start,sORF$End,sORF$CDS_start,sORF$CDS_end) -> ORF_type_reassign
  # 添加回去
  sORF$ORF_type <- ORF_type_reassign
  # 进一步合并，并根据转录本的类型生成子类
  ## 子类1
  sORF %>% mutate(ORF_type_1=case_when(
    ORF_type=="Truncation" | ORF_type=="Extension" | ORF_type=="Readthrough" | ORF_type=="Other_CDS_variant" ~ "CDS_variant",
    TRUE ~ ORF_type)) -> sORF
  ## 子类2
  sORF %>% mutate(ORF_type_2=case_when(
    grepl("pseudogene",sORF$Transcript_type) ~ "ncORF_pseudogene",
    Transcript_type == "lncRNA" | Transcript_type == "retained_intron" | Transcript_type == "processed_transcript" ~ "ncORF_lncRNA",
    TRUE ~ ORF_type)) -> sORF
  sORF %>% mutate(ORF_type_2=case_when(
    ORF_type_2=="ncORF" ~ "ncORF_other",
    TRUE ~ ORF_type_2)) -> sORF
  ## 子类3
  sORF %>% mutate(ORF_type_3=case_when(
    ORF_type_2=="ncORF_lncRNA" ~ Transcript_type,
    TRUE ~ ORF_type_2)) -> sORF
  sORF %>% mutate(ORF_type=NULL) -> sORF
  return(sORF)
}


###### 更新ORF_id_trans #####
update_orf_id_trans <- function(df){
  library(dplyr)
  correct_incorrect_map_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/correct_incorrect_map.rds"
  correct_incorrect_map <- readRDS(correct_incorrect_map_path)
  merge(correct_incorrect_map,df,by.x="ORF_id_trans_incorrect",by.y="ORF_id_trans") -> df_new_orf_id_trans
  df_new_orf_id_trans %>% dplyr::mutate(ORF_id_trans_incorrect=NULL) %>% dplyr::rename(ORF_id_trans=ORF_id_trans_correct) -> df_new_orf_id_trans_1
  return(df_new_orf_id_trans_1)
}

update_orf_id_trans_dt <- function(df) {
  library(data.table)
  
  # 确保输入是data.table
  if(!is.data.table(df)) setDT(df)
  
  # 加载映射表（考虑缓存）
  correct_incorrect_map_path <- "/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/correct_incorrect_map.rds"
  correct_incorrect_map <- readRDS(correct_incorrect_map_path)
  setDT(correct_incorrect_map)
  
  # 使用data.table合并
  result <- correct_incorrect_map[df, 
                                  on = .(ORF_id_trans_incorrect = ORF_id_trans),
                                  nomatch = NULL]
  
  # 重命名列（原地修改避免复制）
  setnames(result, "ORF_id_trans_correct", "ORF_id_trans")
  result[, ORF_id_trans_incorrect := NULL]
  
  return(result)
}

###### 得到特异性谱图数目和肽段数目 #####
get_ms_info <- function(psm){
  psm %>% count(Protein) -> sep_unique_psm_n
  psm %>% distinct(Peptide,.keep_all = T) %>% count(Protein) -> sep_unique_pep_n
  merge(sep_unique_psm_n,sep_unique_pep_n,by="Protein") %>% rename(Unique_psm_n=n.x,
                                                                   Unique_peptide_n=n.y,
                                                                   ORF_id_trans=Protein) -> ms_info
  return(ms_info)
}
