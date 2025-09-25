#!/bin/bash

################################################
# Script Name: process_uniprot_small_proteins.sh
# Author: LiTing
# Created Time: 2025-05-28
# Last Modified: 2025-01-27
# Version: 2.0
################################################

# Description:
#   This script processes human UniProt proteins ≤150 amino acids and matches them
#   with Ensembl database for comprehensive annotation analysis.
#
# Main Functions:
#   1. Extract all UniProt proteins ≤150aa from human proteome
#   2. Remove duplicate protein sequences
#   3. Match proteins with Ensembl database using PeptideMatch
#   4. Classify proteins by type (Isoform vs PE evidence level)
#   5. Generate comprehensive annotation files
#
# Input Files:
#   - UniProt human protein database (defined by define_annotation_gencode_v41_human)
#   - Ensembl GTF annotation file
#   - Ensembl protein database
#
# Output Files:
#   - perfect_match_to_ensId.uniq.txt: Perfect matches between UniProt and Ensembl
#   - uniprot.human.sep.15.gpe: Gene prediction file for matched proteins
#   - uniprot.human.sep.id.type.txt: Protein ID and type classification
#   - processing.log: Detailed execution log
#
# Dependencies:
#   - seqkit: Sequence processing toolkit
#   - PeptideMatchCMD_1.1.jar: Protein sequence matching tool
#   - define_annotation_gencode_v41_human: Environment setup script
#
# Usage:
#   bash process_uniprot_small_proteins.sh
#
# Notes:
#   - Requires GENCODE v41 human annotation environment
#   - Output directory: ../processed/uniprot_id_ens_id_PE_type/
#   - All processing steps are logged to processing.log
################################################
set -eo pipefail

# Configuration
output_path="$(realpath ../processed/uniprot_id_ens_id_PE_type)"
translate_gtf_script="/home/user/data3/lit/project/sORFs/01-ribo-seq/S3.0c.Uni.translate_gtf.v2.20250325.sh"
log_file="${output_path}/processing.log"
tmp_dir="${output_path}/tmp"

# Initialize
mkdir -p "$tmp_dir"
echo "=== Processing started at $(date) ===" | tee "$log_file"

# Function definitions
pep_match() {
    local jar="/home/user/data2/lit/bin/PeptideMatchCMD_1.1.jar"
    local query="$1"
    local db="$2"
    local db_name="$3"
    
    echo "Running peptide matching for $query against $db..." | tee -a "$log_file"
    [ -d "$db_name" ] || java -jar "$jar" -a index -d "$db" -i "$db_name"
    java -jar "$jar" -a query -i "$db_name" -Q "$query" -o out_fasta.txt
}

# Step 1: Extract UniProt proteins ≤150aa
echo "Step 1: Extracting UniProt proteins ≤150aa..." | tee -a "$log_file"
source /home/user/data2/lit/bin/lit_utils.sh
define_annotation_gencode_v41_human
echo "UniProt protein database: $uniprot_prot_fasta" | tee -a "$log_file"

seqkit seq -n -M 150 "$uniprot_prot_fasta" | cut -f 1 -d " " > "$tmp_dir/uniprot.human.sep.id.txt"
total_proteins=$(wc -l < "$tmp_dir/uniprot.human.sep.id.txt")
echo "Total proteins ≤150aa: $total_proteins" | tee -a "$log_file"

# Step 2: Process protein sequences
echo "Step 2: Processing protein sequences..." | tee -a "$log_file"
cd "$tmp_dir"

# Extract sequences
seqkit grep -f uniprot.human.sep.id.txt "$uniprot_prot_fasta" | seqkit seq -w 0 > uniprot.human.sep.fa

# Remove duplicates
seqkit rmdup -s uniprot.human.sep.fa > uniprot.human.sep.rmdup.fa
unique_proteins=$(seqkit stats -T uniprot.human.sep.rmdup.fa | tail -n +2 | cut -f 4)
echo "Unique proteins after deduplication: $unique_proteins" | tee -a "$log_file"

# Create sequence table
seqkit fx2tab uniprot.human.sep.rmdup.fa | awk -v OFS='\t' '{print $1, $NF}' > uniprot.human.sep.rmdup.tab

# Step 3: Match with Ensembl
echo "Step 3: Matching with Ensembl database..." | tee -a "$log_file"
bash "$translate_gtf_script" "$gtf" "$fa"
pep_match uniprot.human.sep.rmdup.fa prot.fa ensembl.prot &> match.info.txt

# Analyze matching results
no_match_count=$(grep -c "No match" out_fasta.txt || echo "0")
perfect_match_count=$(awk '$4==1 && $3==$5' out_fasta.txt | sort -u -k1,1 | wc -l)

echo "Matching statistics:" | tee -a "$log_file"
echo "  - No match: $no_match_count" | tee -a "$log_file"
echo "  - Perfect match: $perfect_match_count" | tee -a "$log_file"
echo "  - Match rate: $(echo "scale=2; $perfect_match_count * 100 / $unique_proteins" | bc)%" | tee -a "$log_file"

# Step 4: Generate output files
echo "Step 4: Generating output files..." | tee -a "$log_file"
awk '$4==1 && $3==$5' out_fasta.txt | sort -u -k1,1 > ../perfect_match_to_ensId.uniq.txt
grep -f <(cut -f2 ../perfect_match_to_ensId.uniq.txt) "$gpe_15" > ../uniprot.human.sep.15.gpe

# Step 5: Classify protein types
echo "Step 5: Classifying protein types..." | tee -a "$log_file"
grep -f "$tmp_dir/uniprot.human.sep.id.txt" "$uniprot_prot_fasta" > "$tmp_dir/uniprot.human.sep.full.id.txt"

awk '
{
    if ($0 ~ /Isoform/) {
        print "Isoform"
    } else {
        match($0, /PE=([0-9]+)/, arr)
        pe_value = (arr[1] != "" ? arr[1] : 1)
        print "PE_" pe_value
    }
}' "$tmp_dir/uniprot.human.sep.full.id.txt" > "$tmp_dir/uniprot.human.sep.type.txt"

paste -d '\t' "$tmp_dir/uniprot.human.sep.id.txt" "$tmp_dir/uniprot.human.sep.type.txt" > ../uniprot.human.sep.id.type.txt

# Step 6: Summary statistics
echo "Step 6: Generating summary statistics..." | tee -a "$log_file"
isoform_count=$(grep -c "Isoform" "$tmp_dir/uniprot.human.sep.type.txt" || echo "0")
pe_count=$(grep -c "PE_" "$tmp_dir/uniprot.human.sep.type.txt" || echo "0")

echo "Protein type classification:" | tee -a "$log_file"
echo "  - Isoform proteins: $isoform_count" | tee -a "$log_file"
echo "  - PE proteins: $pe_count" | tee -a "$log_file"

# Output file summary
echo "Output files generated:" | tee -a "$log_file"
echo "  - ${output_path}/perfect_match_to_ensId.uniq.txt" | tee -a "$log_file"
echo "  - ${output_path}/uniprot.human.sep.15.gpe" | tee -a "$log_file"
echo "  - ${output_path}/uniprot.human.sep.id.type.txt" | tee -a "$log_file"
echo "  - ${output_path}/processing.log" | tee -a "$log_file"

echo "=== Processing completed at $(date) ===" | tee -a "$log_file"