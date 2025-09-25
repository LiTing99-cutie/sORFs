ln -s /home/user/data2/lit/project/ZNF271/01-ribo-seq ./
cp /home/user/data3/rbase/denovo_tumor/denovo_genes/ribosome_profiling/reanalyzed/Chothani_2022_MolCell-brain_embryonic/run.ribo-seq.RiboTISH.sh ./

# 2022_Nature Biotechnology_Mudge et al_Standardized annotation of translated open reading frames
gencode-riboseqORFs

# RiboCode compare
ORFcalling-master

# 2024-BIB
Ribo-seq-Tool-Comparison-Scripts-v2.0