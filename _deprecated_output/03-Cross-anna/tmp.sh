source activate base
# sorf_list=$project_path/02-Mass-spec/human/new_sep_list.txt
sorf_list=$1
output_path=$2
dir=$3
project_path=/home/user/data3/lit/project/sORFs
today=$(date +"%Y%m%d")
mkdir -p $output_path/$dir && cd $output_path/$dir
generate_sorf_genepred_script=$project_path/03-Cross-anna/Uni.gen.genepred.i_sorf_id.py
translate_script=$project_path/01-ribo-seq/S3.0c.Uni.translate_gtf.v2.20250325.sh
source /home/user/data2/lit/bin/lit_utils.sh
define_annotation_gencode_v41_human
python $generate_sorf_genepred_script \
	$sorf_list \
	$gpe_15 \
	gpe
genePredToGtf file gpe gtf
bash $translate_script gtf $fa
seqkit fx2tab cds.fa > cds.tab
seqkit fx2tab prot.fa > prot.tab
echo -e "ORF_id_trans\tSeq\tScodon" > header.txt
paste prot.tab <(less cds.tab |awk '{print substr($2, 1, 3)}') > main.txt
cat header.txt main.txt > sorf.id.seq.scodon.txt