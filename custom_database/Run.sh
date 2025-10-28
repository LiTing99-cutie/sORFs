seqkit fx2tab -l -n -i -H uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta.gz > uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.length.txt

# 验证一下是否正确
seqkit grep -n -r -p "DT3UO_MOUSE" uniprotkb_Mus_musculus_reviewed_canonical_isoform_2024_10_25.fasta.gz > DT3UO_MOUSE.fasta
# MLKMSGWQRQSQNNSRNLRRECSRRKCIFIHHHT