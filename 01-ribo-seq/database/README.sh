# SmProt2; Download in 20241107
wget -r -l1 -H -nd -A "SmProt2*.gz" -e robots=off http://bigdata.ibp.ac.cn/SmProt/download.htm
## After inspect, the mouse data is based on mm10
## should be liftover
cd database/SmProt2
liftOver SmProt2_mouse_Ribo.mm10.bed /home/user/BGM/lit/liftover/mm10/mm10ToMm39.over.chain.gz SmProt2_mouse_Ribo.mm10Tomm39.bed unMapped
# MetamORF; Download in 20241107
wget -r -l1 -H -nd -A "*zip" -e robots=off https://metamorf.hb.univ-amu.fr/downloads

# Openprot
wget -r -l1 -H -nd -A "*zip" -e robots=off https://api.openprot.org/api/2.0/HS/downloads/
https://api.openprot.org/api/2.0/HS/downloads/human-openprot-2_1-refprots+altprots+isoforms-uniprot2022_06_01.tsv.zip

# MicroProteinDB
mkdir MicroProteinDB
wget -r -l1 -H -nd -A "*txt" -e robots=off http://bio-bigdata.hrbmu.edu.cn/MicroProteinDB/download.jsp -P MicroProteinDB
wget -r -l1 -H -nd -A "*fasta" -e robots=off http://bio-bigdata.hrbmu.edu.cn/MicroProteinDB/download.jsp -P MicroProteinDB

# 2022_HongWeiWang_NAR
liftOver nORF.bed /home/user/BGM/lit/liftover/mm10/mm10ToMm39.over.chain.gz nORF.mm10Tomm39.bed unMapped
