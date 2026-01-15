# 0.0 导出de novo ORF所在的基因的gtf，查看目的基因在目的细胞系的RNA-seq数据的表达量
## 需要导出所有的转录本
create_denovo_gtf(){
    source activate biotools
    map=/home/user/data3/lit/project/sORFs/09-CustomDb/formal_20250821/processed/annotation/RibORF_annot/candidate_ORFs/filtered_v2_20251021/representative.tsv
    denovo_list=$1
    output_path=$2
    cut -f1 $denovo_list | tail -n +2 > $output_path/denovo_list.id.txt
    iso_seq_gtf=/home/user/data3/lit/project/sORFs/08-Iso-seq-20250717/results/custom.gtf.with_orf.gtf
    python extract_transcript_gtf.py \
        $output_path/denovo_list.id.txt \
        $map \
        $iso_seq_gtf \
        $output_path/
    # 得到$output_path/denovo_orfs.transcripts.complete.gtf
    grep -w "gene_id" $output_path/denovo_orfs.transcripts.complete.gtf| awk -F 'gene_id "' '{print $2}' | cut -d '"' -f1 | sort -u | wc -l
}

create_denovo_gtf ../processed/new_list_20251128/denovo_list.txt ../processed/new_list_20251128/