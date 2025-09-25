source activate sqanti3
mkdir -p ../results/sqanti/output
cd ../results/sqanti
sqanti3 filter -c sqanti3_config.yaml
sqanti3 rescue -c sqanti3_config.yaml