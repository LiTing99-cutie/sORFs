# mm39
FA=/home/user/data3/lit/resource/genome/mouse/mm39/mm39.fa
GTF=/home/user/data3/lit/resource/gtf/mouse/mm39/gencode.vM29.annotation.gtf
prepare_transcripts -g $GTF -f $FA -o ./

# hg38
source activate ribocode
FA=/home/user/data2/lit/DATA/database/public/genome/hg38/hg38.fa
GTF=/home/user/data2/lit/project/ZNF271/data/annotation/gencode.v41.annotation.gtf
prepare_transcripts -g $GTF -f $FA -o ./