for chr in $(cat ../processed/get_input/chr.list);do
    echo -e "***Processing $chr at $(date '+%Y-%m-%d %H:%M:%S')"    
    bash run.step_2_3.sh $chr
done