################################################
#File Name: run.verify_denovo.sh
#Author: rbase    
#Mail: xiaochunfu@stu.pku.edu.cn
#Created Time: Tue 23 Apr 2024 03:06:51 PM CST
################################################

#!/bin/sh 
#并发运行脚本，并控制并发数
# 设置并发的进程数
thread_num=60
a=$(date +%H%M%S)
# mkfifo
tempfifo="my_temp_fifo"
mkfifo ${tempfifo}
# 使文件描述符为非阻塞式
exec 6<>${tempfifo}
rm -f ${tempfifo}

# 为文件描述符创建占位信息
for ((i=1;i<=${thread_num};i++))
do
{
    echo 
}
done >&6 #事实上就是在fd6中放置了$thread个回车符

workDir=/home/user/data3/rbase/denovo_tumor/denovo_genes/gentree
scriptDir=$workDir/evolution_orfs
gtfFile=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/Homo_sapiens.GRCh38.gencode.v43.annotation.gtf
mafDir=/home/user/data3/rbase/120_mammal_alignment/maf
refDir=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38
refDir=/home/user/data3/rbase/genome_ref/Homo_sapiens/hg38/fasta
genome=$refDir/Homo_sapiens.GRCh38.primary_assembly.genome.fa
transcriptPC=$refDir/transcripts/gencode.v43.pc_transcripts.fa
translationHuman=$refDir/translations/gencode.v43.pc_translations.fa
allPeptidesDir=/home/user/data3/rbase/genome_ref/other_species/peptides
homoDetectDir=$workDir/homolog_detection
outGroupDir=$workDir/outgroup_species
outDir=$workDir/MA_out
ORFDir=$workDir/ORFs_bed
pepDir=$workDir/peptide_fa

if [[ 1 > 2 ]];then
    # find protein-coding de novo genes
    grep -P "\tA" $workDir/hg38_ver95_mechanism.tsv > $workDir/hg38_ver95_mechanism.denovo.tsv
    echo -n > $workDir/denovo.pc.v43.gtf
    cut -f 1 $workDir/hg38_ver95_mechanism.denovo.tsv | while read gene_id;
    do
        grep $gene_id $gtfFile | grep "protein_coding" >> $workDir/denovo.pc.v43.gtf
    done
    # cds gtf
    [ -f $workDir/denovo.pc.v43.cds.gtf ] || grep -P "\tCDS" $workDir/denovo.pc.v43.gtf > $workDir/denovo.pc.v43.cds.gtf
    # cds bed (1-based)
    [ -f $workDir/denovo.pc.v43.cds.bed6 ] || awk 'BEGIN{OFS=FS="\t"}\
            {
                split($9,A,"\"")
                print $1,$4,$5,A[4],".",$7;  
            }' $workDir/denovo.pc.v43.cds.gtf > $workDir/denovo.pc.v43.cds.1-based.bed6
    # cds bed (0-based)
    [ -f $workDir/denovo.pc.v43.cds.bed6 ] || awk 'BEGIN{OFS=FS="\t"}{print $1,$2-1,$3,$4,$5,$6;}' \
        $workDir/denovo.pc.v43.cds.1-based.bed6 > $workDir/denovo.pc.v43.cds.bed6
    # cds fasta
    bedtools getfasta -s -name -fi $refDir/fasta/Homo_sapiens.GRCh38.primary_assembly.genome.fa \
        -bed $workDir/denovo.pc.v43.cds.bed6 \
        -fo $workDir/denovo.pc.v43.cds.fa
    # connect chopped cds fasta sequence
    # 使用awk来处理Fasta文件
    awk '/^>/ { split($0,tA,"::"); starnd=match($0,"(+)"); if (current_id != tA[1]) { current_id = tA[1]; print seq; seq=""; print $0;}} 
        !/^>/ { seq = seq $0 } END { print seq }' $workDir/denovo.pc.v43.cds.fa > $workDir/denovo.pc.v43.cds.whole.fa
    sed -i '1d; s/::chr.*//g' $workDir/denovo.pc.v43.cds.whole.fa
    # remove duplicated CDS (>30nt)
    awk 'BEGIN {RS=">"; ORS=""} NR>1 {sub("\n", "", $2); if (!seen[$2]++ && length($2)>30) print ">" $0}' \
        $workDir/denovo.pc.v43.cds.whole.fa > $workDir/denovo.pc.v43.cds.uniq.fa
    # Translate DNA .fa file to peptide .fa
    faTrans $workDir/denovo.pc.v43.cds.uniq.fa $workDir/denovo.pc.v43.aa.fa
    awk '/^>/ {if (NR!=1) {printf("\n")}; printf("%s\n",$0); next} {printf("%s",$0)} END {printf("\n")}' \
        $workDir/denovo.pc.v43.aa.fa > $workDir/denovo.pc.v43.aa.uniq.fa
    # remove unintact ORF (including Z), retain ATG start codon
    awk 'BEGIN {RS=">"; ORS=""} NR>1 {sub("\n", "", $2); if (index($2, "Z")==0 && match($2, "^M")!=0) print ">" $0}' \
        $workDir/denovo.pc.v43.aa.uniq.fa > $workDir/denovo.pc.v43.aa.intact.fa
    # gene list
    echo "Gene ID" > $workDir/denovo.pc.intact.txt
    grep ">" $workDir/denovo.pc.v43.aa.intact.fa | sed 's/>//g' >> $workDir/denovo.pc.intact.txt

    ######################################################
    ## Get all ORFs (1-based) bed file per de novo gene ##
    ######################################################
    grep ">" $workDir/denovo.pc.v43.aa.intact.fa | sed 's/>//g' | while read tx
    do
        mkdir -p $ORFDir/$tx
        grep $tx $workDir/denovo.pc.v43.cds.1-based.bed6 > $ORFDir/$tx/$tx.ORF.bed
    done

    ######################################################
    ## Get all ORFs peptide fasta file per de novo gene ##
    ######################################################
    grep ">" $workDir/denovo.pc.v43.aa.intact.fa | sed 's/>//g' | while read tx
    do
        mkdir -p $pepDir/$tx
        grep -A 1 $tx $workDir/denovo.pc.v43.aa.intact.fa > $pepDir/$tx/$tx.ORF_pep.fa
    done

    ##################################################################################
    ## Run pipeline to evaluate the evolution of ORF sequences based on 120 mammals ##
    ##################################################################################
    # conda activate bio_seq
    echo "Run pipeline to evaluate the evolution of ORF sequences"
    mkdir -p $scriptDir/tmp/results_120
    mkdir -p $scriptDir/results_120
    mkdir -p $workDir/logs
    species=(hg38 panTro5 gorGor5 ponAbe2 nomLeu3 rheMac8 calJac3 otoGar3)
    cd $scriptDir
    genes=(`tail -n +2 $workDir/denovo.pc.intact.txt`)
    for gene_id in "${genes[@]}";
    do
        echo "### for $gene_id ###"
        # 一个read -u6命令执行一次，就从FD6中减去一个回车符，然后向下执行
        # 当FD6中没有回车符时，就停止，从而实现线程数量控制
        read -u6
        {
            echo "Calculating multiple alignments."
            python3 1_extract_multiple_alignments.py \
                -b $ORFDir/$gene_id/$gene_id.ORF.bed \
                -m $mafDir \
                -o results_120/$gene_id \
                -f yes

            # Calculating ancestral sequences and estimnate intact ORF ancestrally
            echo "Calculating ancestral sequences and estimnate intact ORF ancestrally"
            bash 2.1_ancestral_sequences.sh results_120/$gene_id \
                $workDir/tree/120mammal.nwk \
                $ORFDir/$gene_id/$gene_id.ORF.bed.120mammals \
                $pepDir/$gene_id/$gene_id.ORF_pep.fa

            # Running BLAST across orthologous regions to get conservation score
            # echo "Running BLAST across orthologous regions to get conservation score"
            # python3 3_sequence_conservation_mammals.py -i results_120/$gene_id -p $pepDir/$gene_id.ORF_pep.fa
            
            # Running BLAST across orthologous regions to find specifc proteins
            echo "Running BLAST across orthologous regions to find specifc proteins"
            python3 4_sequence_specific.py --prot_dir results_120/$gene_id --prot_tar $pepDir/$gene_id/$gene_id.ORF_pep.fa

            # select representative primates and re-align
            echo -n > results_120/$gene_id/orfs/$gene_id.dna.fa
            echo -n > results_120/$gene_id/orfs/$gene_id.pep.fa
            for sp in "${species[@]}";
            do 
                grep -A 1 $sp results_120/$gene_id/orfs/$gene_id.maf | sed 's/-//g' >> results_120/$gene_id/orfs/$gene_id.dna.fa
                grep -A 1 $sp results_120/$gene_id/orfs/$gene_id.fa >> results_120/$gene_id/orfs/$gene_id.pep.fa
            done
            mafft --maxiterate 1000 --preservecase --genafpair results_120/$gene_id/orfs/$gene_id.dna.fa \
                > results_120/$gene_id/orfs/$gene_id.closed.maf
            mafft --maxiterate 1000 --preservecase --genafpair results_120/$gene_id/orfs/$gene_id.pep.fa \
                > results_120/$gene_id/orfs/$gene_id.closed.pep.maf
            rm results_120/$gene_id/orfs/$gene_id.dna.fa
            rm results_120/$gene_id/orfs/$gene_id.pep.fa

            # 当进程结束以后，再向FD6中加上一个回车符，即补上了read -u6减去的那个
            echo >&6
        } 1>$workDir/logs/$gene_id.log 2>&1 &
    done
    wait

    ###############################################
    ## Collect evolutionary age of de novo genes ##
    ###############################################
    echo -e "Gene ID\tspeices\tEvolutionary Age\tSyntenic Age\tGained\tGained Convergent\tLost\tDe novo\tSequence" \
        > $workDir/denovo_orf.origination_age.120mammals.txt
    echo -e "Gene ID\tProtein Age\tFarthest Species\tPercent identity\tEvalue\tQuery coverage" \
        > $workDir/denovo_orf.prot_age.txt
    tail -n +2 $workDir/denovo.pc.intact.txt | while read tx
    do
        # collect
        [ -f $ORFDir/$tx/$tx.ORF.bed.120mammals.ancestors ] && \
            tail -n +2 $ORFDir/$tx/$tx.ORF.bed.120mammals.ancestors >> $workDir/denovo_orf.origination_age.120mammals.txt
        [ -f $scriptDir/results_120/$tx/orfs/all.prot_spec.out ] && \
            tail -n +2 $scriptDir/results_120/$tx/orfs/all.prot_spec.out >> $workDir/denovo_orf.prot_age.txt
    done

    #############################
    ## BLAT CDS against genome ##
    #############################
    echo "BLAT CDS against genome"

    # -t=type  Database type.  Type is one of: dna - DNA sequence, prot - protein sequence, dnax - DNA sequence translated in six frames to protein. The default is dna.
    # -q=type  Query type.  Type is one of: dna - DNA sequence, rna - RNA sequence, prot - protein sequence, dnax - DNA sequence translated in six frames to protein, rnax - DNA sequence translated in three frames to protein. The default is dna.
    # -fine  For high-quality mRNAs, look harder for small initial and terminal exons.  Not recommended for ESTs
    # pident > 0.5: -minIdentity=N Sets minimum sequence identity (in percent). Default is 90 for nucleotide searches, 25 for protein or translated protein searches.
    blat $genome $workDir/denovo.pc.v43.cds.uniq.fa \
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
    blastn -query $workDir/denovo.pc.v43.cds.uniq.fa -db $refDir/transcripts/gencode.v43.transcripts \
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
    blastn -query $workDir/denovo.pc.v43.cds.uniq.fa -db $refDir/transcripts/gencode.v43.pc_transcripts \
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

    echo -n > $homoDetectDir/times_in_genome_trans.txt
    genes=(`tail -n +2 $workDir/denovo.pc.intact.txt`)
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

    ##############################################################
    ## BLASTp Peptides against all peptides of outgroup species ##
    ##############################################################
    echo "BLASTp Peptides against all peptides of outgroup species"
    peptides=(`cd $allPeptidesDir && ls *all.fa`)
    for peptide in ${peptides[@]};
    do
        species=`echo $peptide | cut -d "." -f 1`
        echo "---- $species ----"
        # make database for blast
        [ -f $allPeptidesDir/blastp/$species.pdb ] || makeblastdb -in $allPeptidesDir/$peptide -dbtype prot \
            -out $allPeptidesDir/blastp/$species
        # -query sequences.fasta 指定要查询的FASTA文件。
        # -subject 
        # -db nr 指定要搜索的数据库。
        # -out results.txt 指定输出结果的文件名。
        # -outfmt 指定输出格式。上述命令使用格式字符串来指定输出的列。你可以根据需要调整格式字符串
        # -query sequences.fasta 指定要查询的FASTA文件。
        # -subject 
        # -db nr 指定要搜索的数据库。
        # -out results.txt 指定输出结果的文件名。
        # -outfmt 指定输出格式。上述命令使用格式字符串来指定输出的列。你可以根据需要调整格式字符串
        blastp -query $workDir/denovo.pc.v43.aa.intact.fa -db $allPeptidesDir/blastp/$species \
            -out $outGroupDir/denovo_pep.blastp.${species}_peptides.results.txt -num_threads 30 -evalue 1e-5 \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs evalue bitscore" \
            1>blastp.log 2>&1
        # select E-value < 1e-6 && pident > 0.5 && query coverage > 0.5
        awk 'BEGIN{OFS=FS="\t"}\
            {
                if($3>40 && $13>50) print $0;
            }' $outGroupDir/denovo_pep.blastp.${species}_peptides.results.txt > \
            $outGroupDir/denovo_pep.blastp.${species}_peptides.results.real.txt
    done

    ########################################################################################
    ## Count similar proteins in outgroup species of de novo genes and assign origin node ##
    ########################################################################################
    python evolution_orfs/4_origin_outgroup_proteomes.py \
        --gene_file $workDir/denovo.pc.intact.txt \
        --blastp_dir $outGroupDir \
        --output_dir $workDir
fi
    
    ###################################################################
    ## Extract gtf of older de novo genes (not human/hominoid genes) ##
    ###################################################################
    # Rscript /home/user/data3/rbase/denovo_tumor/denovo_genes/gentree/denovo.stats.Rmd
    tx_ids=""
    while read tx_id
    do
        tx_ids+=" -e $tx_id"
    done <<< "$(tail -n +2 $workDir/denovo_genes.older.txt | cut -f 3)"
    grep $tx_ids $workDir/denovo.pc.v43.gtf > $workDir/denovo_genes.older.gtf