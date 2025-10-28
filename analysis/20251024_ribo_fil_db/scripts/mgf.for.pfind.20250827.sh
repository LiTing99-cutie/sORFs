cd /rd1/user/lit/project/sORFs/src
git clone https://github.com/pFindStudio/pUniFind.git
script=/rd1/user/lit/project/sORFs/src/pUniFind/mgf_processor.py
source activate base
find /rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata \
  -name "*mgf" -mtime +30 \
  | grep 21pcw \
  | grep -v -e less3K -e 4_PC \
  | xargs -I{} ln -s {} ../processed/pfind_mgf/

nohup python $script -i ../processed/pfind_mgf/ -o ../processed/pfind_mgf/processed/ --remove-prefixes RTINSECONDS= SCANS= 1/K0= &
