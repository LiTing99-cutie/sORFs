find /rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/CAD20250514licq_BSEP_DDA_60min -name "*Slot2*.d" > file_list.txt
awk -v FS=',' -v OFS='\t' '{print $5,$1}' /rd1/user/lit/project/sORFs/raw_data/MS_metadata_20250523_human.txt > sample_mapping.txt
bash batch_sample_symlink.20250523.sh