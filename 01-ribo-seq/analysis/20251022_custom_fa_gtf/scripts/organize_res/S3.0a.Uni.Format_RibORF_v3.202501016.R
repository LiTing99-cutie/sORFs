# Usage:
# riborf.R <orfs> <gpe> <stat_cutoff> <score_cutoff> <only_ATG> <only_cano_type> <out> [min_aa] [max_aa]

args <- commandArgs(TRUE)
source("/home/user/data2/lit/bin/lit_utils.R"); lib_text()

orfs_path <- args[1]
gpe_path  <- args[2]
stat_cutoff_path <- args[3]
score_cutoff <- args[4]                # custom / youden / fixed
only_canonical_start_codon <- as.numeric(args[5])
only_canonical_type <- as.numeric(args[6])
output_file <- args[7]

# 仅当提供了 min_aa 或 max_aa 时启用长度过滤；否则不启用
has_min <- (length(args) >= 8 && nzchar(args[8]) && args[8] != "NA")
has_max <- (length(args) >= 9 && nzchar(args[9]) && args[9] != "NA")

min_aa <- if (has_min) as.numeric(args[8]) else NA
max_aa <- if (has_max) as.numeric(args[9]) else NA

# 换算为 nt（RibORF 的 length 列为 nt）
min_nt <- if (!is.na(min_aa)) min_aa * 3 + 3 else -Inf
max_nt <- if (!is.na(max_aa)) max_aa * 3 + 3 else  Inf

fread_c(orfs_path) -> orfs
fread_c(gpe_path)  -> genePred

# 拆解 orfID
orfs <- orfs %>%
  tidyr::separate("orfID",
                  into=c("ID","Chr","Strand","Rank","Span","RelaStart","RelaEnd","Type","Scodon"),
                  sep="[|:]", remove=FALSE)

# 长度过滤：仅在提供至少一个边界时启用；否则不筛
if (has_min || has_max) {
  orfs_1 <- subset(orfs, length >= min_nt & length <= max_nt)
} else {
  orfs_1 <- orfs
}

# 类型过滤
if (only_canonical_type) {
  nonCanonical_sorfs <- subset(orfs_1, Type != "canonical")
} else {
  nonCanonical_sorfs <- orfs_1
}

# 起始密码子过滤
if (only_canonical_start_codon) {
  nonCanonical_sorfs <- subset(nonCanonical_sorfs, Scodon == "ATG")
}

# 选择 cutoff
fread_c(stat_cutoff_path) -> A
if (score_cutoff == "custom") {
  B <- dplyr::filter(A, False.pos.rate <= 0.05 & True.pos.rate >= 0.9) |>
       dplyr::filter(True.pos.rate == max(True.pos.rate)) |>
       dplyr::filter(False.pos.rate == min(False.pos.rate))
  cutoff <- B$cutoff
} else if (score_cutoff == "youden") {
  fpr <- A[,6]; tpr <- A[,7]
  cutoff <- A$cutoff[ which.max(tpr - fpr) ]
} else {
  cutoff <- 0.7
}
cat("cutoff is:", cutoff, "\n")

nonCanonical_sorfs <- subset(nonCanonical_sorfs, pred.pvalue >= cutoff)

# 重命名并导出
nonCanonical_sorfs$new_id <- paste0(nonCanonical_sorfs$ID, nonCanonical_sorfs$Strand,
                                     nonCanonical_sorfs$Chr, ":", nonCanonical_sorfs$codon5,
                                     "-", nonCanonical_sorfs$codon3)
m <- merge(nonCanonical_sorfs[, c("orfID","new_id")], genePred, by.x="orfID", by.y="V1")
m <- m[, -1]
fwrite(m, file=output_file, sep='\t', col.names=FALSE)

message(sprintf("[INFO] length_filter=%s; min_aa=%s max_aa=%s -> min_nt=%s max_nt=%s; n=%d",
                ifelse(has_min||has_max, "ON", "OFF"),
                ifelse(has_min, as.character(min_aa), "NA"),
                ifelse(has_max, as.character(max_aa), "NA"),
                ifelse(has_min, as.character(min_nt), "NA"),
                ifelse(has_max, as.character(max_nt), "NA"),
                nrow(m)))
