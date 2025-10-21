output_path=../processed/get_input
mkdir -p $output_path
orfs_gtf=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/sub.gtf
fa=/home/user/data/lit/database/public/genome/hg38/hg38.fa
script=/home/user/data3/lit/project/sORFs/01-ribo-seq/scripts_collapse_20250730/S3.0c.Uni.translate_gtf.v2.20250325.sh
get_cds_bed(){
	less  $1 | awk '$3=="CDS"'| awk -v OFS='\t' '{match($0,/gene_id "([^"]+)"/,arr);print $1,$4-1,$5,arr[1],$6,$7}' > $2
}

orfs_cds_bed=$output_path/orfs.cds.bed
orfs_cds_1_based_bed=$output_path/orfs_cds_1_based.bed

# 得到cds bed
get_cds_bed $orfs_gtf $orfs_cds_bed
awk -v OFS='\t' '{print $1,$2+1,$3,$4,$5,$6}' $orfs_cds_bed > $orfs_cds_1_based_bed
# 得到prot.fa
bash $script $orfs_gtf $fa ../processed/get_input/translate_tmp
orfs_pep_fa=../processed/get_input/translate_tmp/prot.fa
# 得到分染色体的cds bed和prot.fa
chr_list_file="${output_path}/chr.list"
cut -f1 "${orfs_cds_1_based_bed}" | sed 's/^chr/chr/' | sort -V | uniq > "${chr_list_file}"
id_list=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/noncano.id.lst
mkdir -p ${output_path}/per_chr/orf_list/
while read -r chr; do
    grep "\:$chr\:" $id_list > "${output_path}/per_chr/orf_list/${chr}.orf.list"
done < "${chr_list_file}"

while read -r chr; do
    mkdir -p ${output_path}/per_chr/ORFs_bed/$chr ${output_path}/per_chr/peptide_fa/$chr
    for id in $(cat ${output_path}/per_chr/orf_list/${chr}.orf.list);do
        grep -F "$id" $orfs_cds_1_based_bed > ${output_path}/per_chr/ORFs_bed/$chr/$id.ORF.bed
        grep -A1 -F "$id" $orfs_pep_fa > ${output_path}/per_chr/peptide_fa/$chr/$id.ORF_pep.fa
    done 
done < "${chr_list_file}"