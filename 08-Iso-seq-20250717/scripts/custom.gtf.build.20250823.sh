gffread ../processed/collapsed.sorted.filtered_lite.gff -T -o ../results/custom.gtf
gtfToGenePred -geneNameAsName2 -genePredExt ../results/custom.gtf ../results/custom.15.gpe