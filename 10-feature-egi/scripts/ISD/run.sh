fa_dir=/home/user/data3/lit/project/sORFs/10-feature-egi/processed/seve_list_compare
aminoAcidFastaPC=$fa_dir/canonical_orfs.sample.2.5k.lt.150aa.fa
aminoAcidFastalncRNA=$fa_dir/lnc_orfs.sample.2.5k.lt.150aa.fa
aminoAcidFastaIntergenic=$fa_dir/intergenic_orfs.sample.2.5k.lt.150aa.aa.fa
aminoAcidFastaNoncano=$fa_dir/non_cano_orfs.sample.2.5k.lt.150aa.fa
aminoAcidFastaNoncano_r1_m_minus=$fa_dir/r1_m_minus_non_cano_orfs.sample.2.5k.lt.150aa.fa

nohup bash ISD/compare_ISD.20251116.sh $aminoAcidFastaPC $aminoAcidFastalncRNA $aminoAcidFastaIntergenic $aminoAcidFastaNoncano $aminoAcidFastaNoncano_r1_m_minus &> ../log/ISD.20251116.log &
nohup bash ISD/compare_ISD.20251116.sh /home/user/data3/lit/project/sORFs/11-denovo-list/processed/new_list_20251128/denovo_list.fa &> ../log/ISD.20251128.log &

cat ../processed/ISD/ISD_avg/*ISD_avg*txt > ../processed/ISD/ISD.merged.txt