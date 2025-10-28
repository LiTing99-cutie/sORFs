stat_path_1 <- "output/stat/stat.rds"
stat_path_2 <- "output/stat/default_trans/stat.rds"
stat_path_3 <- "output/stat/open_ribo/stat.rds"
sample_order_path <- "output/sample_order.rds"
output_path_1 <- "output/stat/total/"

readRDS(stat_path) %>% .$all_sample_all_level_sep %>% 
  filter(Confidence_Level==1) -> df_level_1
readRDS(sample_order_path) -> sample_order
