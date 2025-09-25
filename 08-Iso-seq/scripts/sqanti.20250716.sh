source activate sqanti3
mkdir -p ../results/sqanti/output
cd ../results/sqanti
[ -f collapsed.sorted.gtf ] || gffread ../../processed/collapsed.sorted.gff -T -o collapsed.sorted.gtf
sqanti3 all -c sqanti3_config.yaml