database=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/group_specific_db/uniprot.contam.sorfs.trans.modiHeader.fa
cd output/S1
# 和我们鉴定出的SEP做比较，看看其他人鉴定出的peptide是否支持我们的SEP
# 得到我们的SEP list
less /home/user/data3/lit/project/sORFs/03-Cross-anna/analysis/annotate_ms_orfs_20250401/output/sep_add_basic_ms_ribo_info_group_retained.txt | \
 cut -f 1,4|tail -n +2|seqkit tab2fx > in_house_sep.fa
# 去掉N端的M
sed '{ s/^M// }' in_house_sep.fa > in_house_sep.clipM.fa
# 分别进行酶切消化，去重后合并消化得到的肽段
conda activate base
mkdir tmp
python3 ../../1.1a.uni.peptide.digest.py -i in_house_sep.fa -o tmp/in_house_sep.digested.peptide.fa
python3 ../../1.1a.uni.peptide.digest.py -i in_house_sep.clipM.fa -o tmp/in_house_sep.clipM.digested.peptide.fa
cat tmp/in_house_sep.digested.peptide.fa tmp/in_house_sep.clipM.digested.peptide.fa |seqkit rmdup -s > in_house_sep.digested.peptide.fa
# 看哪些肽段可以完美比对到SEP消化得到的肽段
pep_match(){
    # 输入，每行一个肽段
    # I/L equivalent
        # -a,--action        The action to perform ("index" or "query").
    # -d,--dataFile      The path to a FASTA file to be indexed.
    # -e,--LeqI               Treat Leucine (L) and Isoleucine (I) as
    #                         equivalent (default: no).
    # -f,--force              Overwrite the indexDir (default: no).
    # -h,--help               Print this message.
    # -i,--indexDir      The directory where the index is stored.
    # -l,--list               The query peptide sequence file is a list of
    #                         peptide sequences, one sequence per line
    #                         (default: no).
    # -o,--outputFile    The path to the query result file.
    # -Q,--queryFile     The path to the query peptide sequence file in
    #                         either FASTA format or a list of peptide
    #                         sequences, one sequence per line.
    # -q,--query         One peptide sequence or a comma-separated list of
    jar=/home/user/data2/lit/bin/PeptideMatchCMD_1.1.jar
    query=$1
    db=$2
    db_name=$3
    output=$4
    [ -d $db_name ] || java -jar $jar -a index -d $db -i $db_name
    java -jar $jar -a query -i $db_name -Q $query -l -e -o $output
}
pep_match pub_non_HLA_unique.peptide.fa in_house_sep.digested.peptide.fa in_house_sep pub_in_house_map.txt
less pub_in_house_map.txt|grep -v "No match"|awk '$3==$5&&$4==1' > pub_in_house_map_perfect.txt
less pub_in_house_map_perfect.txt | cut -f 1|uniq > pub_in_house_map_perfect.peptide.txt
# 查看这些完美比对到SEP的肽段是否可以比对到其他的蛋白质上面
## 但是没有考虑到漏切
out=pub_in_house_map_perfect.peptide.whetherUniq.txt
[ -f $out ] || rm -rf $out 
for pattern in $(cat pub_in_house_map_perfect.peptide.txt);do
n=$(grep -E  "^M${pattern}|^${pattern}|[RK]${pattern}" "$database"|sort|uniq|wc -l)
echo -e "$pattern\t$n"
done >> $out 

# 消化全库蛋白
path=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Run_for_human_20250227/human_trans_database_20250324/output/group_specific_db/
seqkit split -p 10 $database  # 分割为10个子文件
parallel -j 10 \
  python ../../1.1a.uni.peptide.digest.py -i {} -o digested_parts/{/.}.digested.fa \
  ::: $path/uniprot.contam.sorfs.trans.modiHeader.fa.split/uniprot.contam.sorfs.trans.modiHeader.part_*.fa
cat digested_parts/uniprot.contam.sorfs.trans.modiHeader.part_*.digested.fa > digested_parts/uniprot.contam.sorfs.trans.modiHeader.total.digested.fa

# clip掉N端的M之后再进行消化
sed '{ s/^M// }' $database > database.clipM.fa
seqkit split -p 10 database.clipM.fa
parallel -j 10 \
  python ../../1.1a.uni.peptide.digest.py -i {} -o digested_parts/{/.}.digested.fa \
  ::: database.clipM.fa.split/database.clipM.part_*.fa
cat digested_parts/database.clipM.part_*.digested.fa > digested_parts/database.clipM.total.digested.fa

# 合并结果，以蛋白质分组，得到去重之后的peptide
cat digested_parts/database.clipM.total.digested.fa  digested_parts/uniprot.contam.sorfs.trans.modiHeader.total.digested.fa > digested_parts/total.fa
rm -rf digested_parts/*part*
seqkit fx2tab digested_parts/total.fa > total.tab
less total.tab | awk '{print $1,$NF}'| \
awk '
{
    # 以第一列作为分组键，第二列作为值
    if (!seen[$1, $2]++) {
        print $0  # 如果该组合尚未出现过，则打印整行
    }
}
' > total.dedup.tab
awk -v OFS='\t' '{print $1,$2}' total.dedup.tab > total.dedup.tab.1 && rm -rf total.dedup.tab
seqkit tab2fx total.dedup.tab.1 >  total.dedup.fa

pep_match pub_non_HLA_unique.peptide.fa total.dedup.fa in_house_db pub_in_house_db_map.txt

less pub_in_house_db_map.txt|grep -v "No match"|awk '$3==$5&&$4==1' > pub_in_house_db_map_perfect.txt
awk '{
    arr[$1][$2]
} 
END {
    for (key in arr) {
        count = 0
        for (id in arr[key]) { count++ }
        print key, count
    }
}' pub_in_house_db_map_perfect.txt > pub_in_house_db_map_perfect_map_n.txt

