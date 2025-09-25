
source "/home/user/data2/lit/bin/lit_utils.sh"
define_annotation_gencode_v41_human
parse_offdict_script=/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250714/parse_offdict.py
script_dir="/home/user/data3/lit/project/sORFs/01-ribo-seq/analysis"
source activate ribocode
extract_specific_length_reads_script="$script_dir/Test-20250509/extract_specific_length_reads.sh"
for bam in $(find "$PWD" -name "*Aligned.sortedByCoord.out.bam"); do
    sample=$(basename $bam | sed 's/_Aligned.sortedByCoord.out.bam//')
    [ -f "$bam.bai" ] || samtools index "$bam"
    outdir=$(dirname "$bam")/ribotish_0.6_cutoff
    mkdir -p $outdir
    # ribotish quality -b "$bam" -g "$gtf" -p 30 --th 0.6 \
    #     -r "$outdir/0.6_qual.txt" \
    #     -f "$outdir/0.6_qual.pdf" \
    #     -o "$outdir/0.6.para.py"
    cd $outdir
    python $parse_offdict_script 0.6_qual.txt > offset.tab.txt
    bash "$extract_specific_length_reads_script" "$bam" "offset.tab.txt" "ribotish.peri.bam"
done

find "$PWD" -name "ribotish.peri.bam" | xargs -I {} sh -c 'cnt=$(samtools view -c "{}"); echo -e "{}\t$cnt"' > /home/user/data3/lit/project/sORFs/01-ribo-seq/analysis/Test-20250714/01-output/stat/p_sites_number_ribotish.0.6.txt

awk -F'/' '{for(i=1;i<=NF;i++) if($i ~ /^p21/) {name=$i}} {print name "\t" $NF}' 01-output/stat/p_sites_number_ribotish.0.6.txt
