denovo_dir=/home/user/data3/lit/project/sORFs/05-denovo-status/analysis/20251029_denovo_check/processed
orf=PB.29415.11__chr22__-__201__3603__2541__2736__dORF__CTG
pep=$denovo_dir/check_denovo/chr22/orfs/$orf.ortholog.pep.fa
cds=$denovo_dir/check_denovo/chr22/orfs/$orf.ortholog.nucl.fa
# -h, -homolog	Homolog group file [string, required]
# -n, -nuc	File containing multiple nucleotide sequences [string, required]
# -a, -aminoacid	File containing multiple amino acid sequences [string, required]
# -m, -msa	Multiple sequence aligner or specify with its full path [clustalw2 | t_coffee | mafft | muscle, optional]
# -f, -format	Output file format [axt | fasta | paml | codon | clustal, optional], default = fasta
# -g, -nogap	Remove aligned codons with gaps
# -t, -nomismatch	Remove mismatched codons
# -k, -kaks	Enable using KaKs_Calculator for Ka and Ks estimation (requiring axt format)
# -o, -output	Output folder [string, required]
out_dir=$(realpath ../processed/evolution_para/kaks/$orf)
mkdir -p $out_dir/paraAT_res
seqkit seq -n $pep |tr "\n" "\t" > $out_dir/homolog
# 生成所有组合
seqkit seq -n $pep | awk '{ids[NR]=$0} END {
    for(i=1; i<=NR-1; i++) {
        for(j=i+1; j<=NR; j++) {
            print ids[i]"\t"ids[j]
        }
    }
}' > $out_dir/homolog
# 人类 vs 其他物种
seqkit seq -n $pep | awk 'NR==1{hg38=$0; next} {print hg38"\t"$0}' > $out_dir/homolog

proc_n=10
echo $proc_n > proc
[ -d $out_dir/paraAT_res ] && rm -rf $out_dir/paraAT_res
ParaAT.pl -h $out_dir/homolog -n $cds -a $pep -processor $PWD/proc -m mafft -f axt -g -t -k -o $out_dir/paraAT_res
# 	-m	Methods for estimating Ka and Ks and theirs references(Default = MA)
		#   NG		Nei, M. and Gojobori, T. (1986) Mol. Biol. Evol., 3, 418-426.
KaKs_Calculator -i $out_dir/*.axt -o $out_dir/*.axt.ng -m NG

echo "=== 计算Ka/Ks ==="
for axt in $out_dir/paraAT_res/*.axt; do
    base=$(basename $axt .axt)
    KaKs_Calculator -i $axt -o $out_dir/paraAT_res/${base}.kaks -m NG
done

echo "=== Ka/Ks结果 ==="
echo -e "Species_pair\tKa\tKs\tKa/Ks\tS-Substitutions\tN-Substitutions\tSelection" >> $out_dir/kaks.res
for kaks in $out_dir/paraAT_res/*aln.kaks; do
    tail -n +2 $kaks | awk '{
        type = ($5 > 1) ? "Positive" : ($5 < 1) ? "Purifying" : "Balancing"
        print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$12"\t"$13"\t"type
    }'
done >> $out_dir/kaks.res