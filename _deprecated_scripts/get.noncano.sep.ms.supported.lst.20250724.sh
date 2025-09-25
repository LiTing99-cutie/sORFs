psm_sep_all=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250723_analysis/MS_res_from_Galaxy/psm_sep_all.txt
output_path=../processed/tab_info
# 得到非经典小肽的ID
cat $psm_sep_all|awk -v FS='\t' -v OFS='\t' '{print $35}'|grep -v sp|sort -u|grep -v Protein > $output_path/noncano.sep.ms.specific.nonspecific.tmp.lst
# 添加一列列名，并删除临时文件
sed '1i ORF_id_trans' $output_path/noncano.sep.ms.specific.nonspecific.tmp.lst > $output_path/noncano.sep.ms.specific.nonspecific.lst
rm -rf $output_path/noncano.sep.ms.specific.nonspecific.tmp.lst
