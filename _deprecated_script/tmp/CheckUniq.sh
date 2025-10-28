#!/usr/bin/sh

################################################
#File Name: CheckUniq.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: 2024年11月20日 星期三 15时06分01秒
################################################

set -eo pipefail

psm=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024_11_14_batch_run/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/psm.tsv

cd debug
less $psm | awk -F'\t' '$33=="true"' | awk -F'\t' '{printf "%s\n%s\n", ">"$34":"$3, $3}' > unique.psm.p.txt
seqkit rmdup -n unique.psm.p.txt > unique.psm.p.rmdup.txt

fa=/rd1/user/lit/project/sORFs/custom_database/uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.nonCano.sorf.filtered.fa
makeblastdb -in $fa -dbtype prot -out ./database

blastp -query unique.psm.p.rmdup.txt -db ./database -outfmt '6 qseqid sseqid pident qlen slen length bitscore evalue' -num_threads 30 -out unique.psm.p.blastp.txt

grep -A 1 "sp|P23116|EIF3A_MOUSE:RGMDDDRGPR" unique.psm.p.rmdup.txt
grep -A 1 "sp|P23116|EIF3A_MOUSE" $fa

cd tmp_1
psm=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024-11-19/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562/psm.tsv
less $psm | awk -F'\t' '$31=="true"' | awk -F'\t' '{printf "%s\n%s\n", ">"$32":"$3, $3}' > unique.psm.p.txt
seqkit rmdup -n unique.psm.p.txt > unique.psm.p.rmdup.txt
fas=/rd1/user/lit/project/sORFs/custom_database/2024-10-31-decoys-reviewed-isoforms-contam-UP000000589.fas
seqkit grep -r -p "^sp" $fas > /rd1/user/lit/project/sORFs/custom_database/UP000000589.fa
fa=/rd1/user/lit/project/sORFs/custom_database/UP000000589.fa
makeblastdb -in $fa -dbtype prot -out ./database
blastp -query unique.psm.p.rmdup.txt -db ./database -outfmt '6 qseqid sseqid pident qlen slen length' -num_threads 30 -out unique.psm.p.blastp.txt

rep_1_ion=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/F2_rerun_new_para/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/ion.tsv
rep_2_ion=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/2024_11_14_batch_run/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/ion.tsv
rep_1_protein=/rd1/user/lit/project/sORFs/output/MS/Fragpipe_output/F2_rerun_new_para/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/protein.tsv
rep_2_protein=/rd1/user/lit/project/sORqFs/output/MS/Fragpipe_output/2024_11_14_batch_run/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_F2_Slot2_52_1_3562_d/protein.tsv

grep AAAASAAEAGIATPGTEDSDDALLK $rep_1_ion
grep AAAASAAEAGIATPGTEDSDDALLK $rep_2_ion

grep O35226 $rep_1_protein
grep O35226 $rep_2_protein

test_path=./output/MS/Fragpipe_output/2024_11_18_batch_run_6_samples_test_overlap/Trypsin_LysC/CAD20241022licq_BSEP_PreExp_DDA_60min_S4_Slot1_54_1_3559_d/
grep P27870 $test_path/protein.tsv
grep P27870 $test_path/ion.tsv
grep P27870 $test_path/peptide.tsv
grep P27870 $test_path/psm.tsv