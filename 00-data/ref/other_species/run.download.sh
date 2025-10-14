# workDir=/home/user/data3/lit/project/sORFs/00-data/ref/other_species
# cat $workDir/cdna.list.txt | while read url;
# do
#     mkdir -p $workDir/cdna/ && cd $workDir/cdna/
#     wget -c $url
# done

# cat $workDir/dna.list.txt | while read url;
# do
#     mkdir -p $workDir/dna/ && cd $workDir/dna/
#     wget -c $url
# done

#!/usr/bin/env bash
set -euo pipefail

workDir="/home/user/data3/lit/project/sORFs/00-data/ref/other_species"
cdna_list="$workDir/cdna.list.txt"
dna_list="$workDir/dna.list.txt"

cdna_dir="$workDir/cdna"
dna_dir="$workDir/dna"

# 并行度：按带宽/CPU调整
P=8

mkdir -p "$cdna_dir" "$dna_dir"

# cDNA
grep -vE '^\s*(#|$)' "$cdna_list" \
  | xargs -n1 -P "$P" -I{} wget -c -nv -P "$cdna_dir" "{}"

# DNA
grep -vE '^\s*(#|$)' "$dna_list" \
  | xargs -n1 -P "$P" -I{} wget -c -nv -P "$dna_dir" "{}"
