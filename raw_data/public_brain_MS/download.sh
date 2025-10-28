#!/usr/bin/sh

################################################
#File Name: download.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2025年04月30日 星期三 16时54分23秒
################################################

set -eo pipefail

less HLA/exp_contrib_table_guest_20250430-014349-995.tsv

tr -d '\r' < non-HLA/exp_contrib_table_guest_20250430-012600-968.tsv > non-HLA/exp_contrib_table_guest_20250430-012600-968.1.tsv
tr -d '\r' < HLA/exp_contrib_table_guest_20250430-014349-995.tsv > HLA/exp_contrib_table_guest_20250430-014349-995.1.tsv

less HLA/exp_contrib_table_guest_20250430-014349-995.1.tsv|cut -f 20| sort|uniq -c
less non-HLA/exp_contrib_table_guest_20250430-012600-968.1.tsv|cut -f 21| sort|uniq -c

less HLA/exp_contrib_table_guest_20250430-014349-995.1.tsv|awk  -F'\t' '$20=="Brain"'|cut -f 1|sort|uniq > HLA/brain_pxd.lst
less non-HLA/exp_contrib_table_guest_20250430-012600-968.1.tsv|awk  -F'\t' '$21=="Brain"'|cut -f 1|sort|uniq > non-HLA/brain_pxd.lst

wget -r -nH --cut-dirs=5 -np --reject="index.html*" -P non-HLA ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD035950/
# wget -r -nH --cut-dirs=5 -np --reject="index.html*" -P HLA ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/
wget -r -nH --cut-dirs=5 -np --reject="index.html*" --accept="*Brain*.raw,*Cerebellum*.raw" -P HLA ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/
# 需要指定下载日期，而不是直接使用PXD
# parallel -j 8 "wget -r -nH --cut-dirs=5 -np --reject="index.html*" -P non-HLA ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/{}/" :::: non-HLA/brain_pxd.lst

# 直接wget从对应PXD上面下载会存在下载不全的问题，需要手动check
cd HLA/PXD019643/
less README.txt |grep Cerebellum.*raw |cut -f 2 > Cerebellum.raw.full.lst
sed 's/%//' Cerebellum.raw.full.lst > Cerebellum.raw.full.1.lst
ls |grep Cerebellum.*raw$ > Cerebellum.raw.download.1.lst
 diff <(sort -k1,1 Cerebellum.raw.download.1.lst) <(sort -k1,2 Cerebellum.raw.full.1.lst )
grep 170107_AM_AUT01-DN13_Cerebellum_W6-32_10%_DDA_1_400-650mz_msms13_standard.raw README.txt
wget ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170107_AM_AUT01-DN13_Cerebellum_W6-32_10_DDA_1_400-650mz_msms13_standard.raw

ls |grep Cerebellum.*raw$ > Cerebellum.raw.download.2.lst
 diff <(sort -k1,1 Cerebellum.raw.download.2.lst) <(sort -k1,2 Cerebellum.raw.full.1.lst )

less README.txt |grep Brain.*raw |cut -f 2 > Brain.raw.full.lst
sed 's/%//' Brain.raw.full.lst > Brain.raw.full.1.lst
ls |grep Brain.*raw$ > Brain.raw.download.1.lst
 diff <(sort -k1,1 Brain.raw.download.1.lst) <(sort -k1,2 Brain.raw.full.1.lst )


wget -r -nH --cut-dirs=5 -np --reject="index.html*" --accept="*Brain*.raw,*Frontalcortex*.raw" -P non-HLA ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/04/PXD000561
# 直接wget从对应PXD上面下载会存在下载不全的问题，需要手动check
ls *raw > download.lst
cut -f 2 README.txt  |grep raw|grep -E "Brain|Frontalcortex" > real.lst
diff <(sort -k1,1 download.lst ) <(sort -k1,1 real.lst )
# 下次可以下载了README.txt之后从中获取链接再进行下载