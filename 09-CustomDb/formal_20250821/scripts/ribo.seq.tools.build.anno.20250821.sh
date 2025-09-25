fa=/home/user/data3/lit/project/sORFs/07-Genome/results/custom_fa/custom_ref.fa
gtf=/home/user/data3/lit/project/sORFs/09-CustomDb/Test_20250801/processed/mkAnno_for_moPepGen/custom.gtf
genePred=/home/user/data3/lit/project/sORFs/09-CustomDb/Test_20250801/processed/mkAnno_for_moPepGen/custom.15.gpe
RibORF_script_path=/home/user/data2/lit/software/RibORF/RibORF.2.0
out_path=../processed/annotation
RiboCode_annot_path=$out_path/RiboCode_annot
RibORF_annot_path=$out_path/RibORF_annot
PRICE_anno_name=hg38_custom

mkdir -p $RiboCode_annot_path
mkdir -p $RibORF_annot_path

# ribocode
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start ribocode"
source activate ribocode
prepare_transcripts -g $gtf -f $fa -o $RiboCode_annot_path
# riborf
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start riborf"
perl $RibORF_script_path/ORFannotate.pl -g $fa -t $genePred -s ATG/CTG/GTG/TTG/ACG -l 21 -o $RibORF_annot_path
# price
echo "[`date '+%Y-%m-%d %H:%M:%S'`] Start price"
gedi -e IndexGenome -s $fa -a $gtf -n $PRICE_anno_name
