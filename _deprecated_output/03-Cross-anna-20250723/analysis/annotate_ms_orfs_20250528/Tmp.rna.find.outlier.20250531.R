
fread_c("./output/S2/rna_seq_qc/Uniquely_mapped_reads_rate_number.txt") -> rna_uniq_n_r
SRR_lst <- c("SRR15906371","SRR15906388","SRR15906382","SRR15906506","SRR15906495",
             "SRR15906389","SRR15906387","SRR15906391","SRR15906372","SRR15906377")
filter(rna_uniq_n_r,Sample %in% SRR_lst)
