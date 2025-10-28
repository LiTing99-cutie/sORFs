```{r}
source("~/bin/lit_utils.R")
source("/rd1/user/lit/project/sORFs/sORFs.utils.R")
read_psm("/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/Formal/mouse/2025_02_09_default_trans_database/Trypsin/CAD20241224licq_BSEP_DDA_60min_E16_1_PAGE_T_T_Slot2_49_1_4736_d/psm.tsv") -> psm
colnames(psm)
```

```{r}
head(psm)
filter(psm,`Is.Unique`=="true") -> unique_psm
unique_psm %>% filter(grepl("sp",Protein)) -> anno_pro_unique_psm
unique_psm %>% filter(grepl("ENS",Protein)) -> sorf_unique_psm

cols <- c("Delta.Mass","RTScore","Ion.Mobility","Hyperscore")
anno_pro_unique_psm[,cols] %>% colMeans()
sorf_unique_psm[,cols] %>% colMeans()
```
```{r}
sample_n <- 20
set.seed(123)
sample(sorf_unique_psm$Spectrum,size = sample_n) %>% cat(sep = "\n")
```

```{r}
read.
```


