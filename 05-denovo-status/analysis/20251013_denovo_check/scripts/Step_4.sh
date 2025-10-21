source activate base
output_path=../processed/blastp_check
mkdir -p $output_path && cd $output_path

# 准备pep以及cds的序列
total_pep_fa=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/translate_out/prot.fa
total_cds_fa=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/augment_orf_table/translate_out/cds.fa
gene_lst=/home/user/data3/lit/project/sORFs/02-Mass-spec-20250723/analysis/20250910_c8_protein_map/results/1/noncano.id.lst
seqkit grep -f $gene_lst $total_pep_fa > target.pep.fa
seqkit grep -f $gene_lst $total_cds_fa > target.cds.fa

# 首先对外类群进行blastp
query=target.pep.fa
db_path=/home/user/data3/rbase/genome_ref/other_species/peptides/blastp/

for file in $(ls $db_path/*);do
basename $file | cut -f 1 -d "."
done | sort|uniq > out_group.blastp.lst

for species in $(cat out_group.blastp.lst);do
blastp -query $query -db $db_path/$species -out $species.res.out -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs evalue bitscore' -num_threads 20
less $species.res.out | awk '$3>40 && $13>50 && $14 < 1e-5' > denovo_pep.blastp.${species}_peptides.results.real.txt
done

# 其次对人这个物种进行blastn
query=target.cds.fa
genome=/home/user/data/lit/database/public/genome/hg38/hg38.fa
# # ensembl 106
# cdna_fa=/home/user/data/lit/database/public/annotation/gene_and_gene_predictions/Homo_sapiens.GRCh38.cdna.all.fa
# makeblastdb -in $cdna_fa -dbtype nucl
refDir=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/fasta
# transcriptPC=$refDir/transcripts/gencode.v43.pc_transcripts.fa
homoDetectDir=$PWD

#############################
## BLAT CDS against genome ##
#############################
echo "BLAT CDS against genome"

# -t=type  Database type.  Type is one of: dna - DNA sequence, prot - protein sequence, dnax - DNA sequence translated in six frames to protein. The default is dna.
# -q=type  Query type.  Type is one of: dna - DNA sequence, rna - RNA sequence, prot - protein sequence, dnax - DNA sequence translated in six frames to protein, rnax - DNA sequence translated in three frames to protein. The default is dna.
# -fine  For high-quality mRNAs, look harder for small initial and terminal exons.  Not recommended for ESTs
# pident > 0.5: -minIdentity=N Sets minimum sequence identity (in percent). Default is 90 for nucleotide searches, 25 for protein or translated protein searches.
blat $genome $query \
    -t=dna -q=dna -fine -out=psl -minIdentity=50 \
    $homoDetectDir/denovo_cds.blat_genome.results.txt 1>blat.log 2>&1

# select query coverage in query seq > 0.5 && pident > 0.5
tail -n +6 $homoDetectDir/denovo_cds.blat_genome.results.txt | awk 'BEGIN{OFS=FS="\t"}\
    {
        if(($13-$12)/$11>0.5) print $0;
    }' > $homoDetectDir/denovo_cds.blat_genome.results.real.txt

######################################
## BLASTn CDS against transcription ##
######################################
echo "BLASTn CDS against transcription"

# -query sequences.fasta 指定要查询的FASTA文件。
# -subject 
# -db nr 指定要搜索的数据库。
# -out results.txt 指定输出结果的文件名。
# -outfmt 指定输出格式。上述命令使用格式字符串来指定输出的列。你可以根据需要调整格式字符串
blastn -query $query -db $refDir/transcripts/gencode.v43.transcripts \
    -out $homoDetectDir/denovo_cds.blastn.results.txt -num_threads 30 \
    -perc_identity 50 -evalue 1e-6 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs evalue bitscore" 1>blastn.log 2>&1
# select E-value < 1e-6 && pident > 0.5 && query coverage > 0.5
awk 'BEGIN{OFS=FS="\t"}\
    {
        if($13>50) print $0;
    }' $homoDetectDir/denovo_cds.blastn.results.txt > $homoDetectDir/denovo_cds.blastn.results.real.txt

##################################################
## BLASTn CDS against all human pc transcripts  ##
##################################################
echo "BLASTn CDS against all human pc transcripts"
# make database for blast
# makeblastdb -in $transcriptPC -dbtype nucl -out $refDir/transcripts/gencode.v43.pc_transcripts
# -query sequences.fasta 指定要查询的FASTA文件。
# -subject 
# -db nr 指定要搜索的数据库。
# -out results.txt 指定输出结果的文件名。
# -outfmt 指定输出格式。上述命令使用格式字符串来指定输出的列。你可以根据需要调整格式字符串
blastn -query $query -db $refDir/transcripts/gencode.v43.pc_transcripts \
    -out $homoDetectDir/denovo_cds.blastn.pc.results.txt -num_threads 30 \
    -perc_identity 50 -evalue 1e-6 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs evalue bitscore" 1>blastn.pc.log 2>&1
# select E-value < 1e-6 && pident > 0.5 && query coverage > 0.5
awk 'BEGIN{OFS=FS="\t"}\
    {
        if($13>50) print $0;
    }' $homoDetectDir/denovo_cds.blastn.pc.results.txt > $homoDetectDir/denovo_cds.blastn.pc.results.real.txt

##########################################################################################################
## Exclude peptides that were mapped to more than 1 genomic locations or transcripts of different genes ##
##########################################################################################################
seqkit seq -n target.pep.fa > target.id.txt
echo -n > $homoDetectDir/times_in_genome_trans.txt
genes=(`cat target.id.txt`)
for gene in ${genes[@]};
do
    echo "--- $gene ---"
    # find if multi-mapping in different locations in blat results
    loc_times_blat=`grep -c $gene $homoDetectDir/denovo_cds.blat_genome.results.real.txt`
    # find if multi-mapping in different genes in blastn results
    gene_times_blastn=`grep $gene $homoDetectDir/denovo_cds.blastn.results.real.txt | cut -f2 | cut -f2 -d"|" | sort -u | wc -l`
    # find if multi-mapping in different genes in blastp results
    gene_times_blastn_pc=`grep $gene $homoDetectDir/denovo_cds.blastn.pc.results.real.txt | cut -f2 | cut -f2 -d"|" | sort -u | wc -l`
    echo "$gene $loc_times_blat $gene_times_blastn $gene_times_blastn_pc" | sed 's/ /\t/g' >> $homoDetectDir/times_in_genome_trans.txt
done

sed '1i Gene ID' target.id.txt > target.id.h.txt
python /home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/evolution_orfs/4_origin_outgroup_proteomes.py \
    --gene_file target.id.h.txt \
    --blastp_dir $PWD \
    --output_dir $PWD