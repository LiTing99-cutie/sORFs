# sample=21pcw_1_C8_T_T
closed_dir=../processed/pFind_res_20250829/closed
open_dir=../processed/pFind_res_20250829/open
meta_file=/rd1/user/lit/project/sORFs/raw_data/In_House_Guomics_SP_MSdata/human_organized_20250625/metadata/metadata.formal.v1.txt
find $closed_dir -name "pFind-Filtered.spectra"|head -n1|xargs cat|head -n1 > $closed_dir/closed.res.header.txt
find $open_dir -name "pFind-Filtered.spectra"|head -n1|xargs cat|head -n1 > $open_dir/open.res.header.txt
for sample in $(cut -f 1 $meta_file -d ','|sed 's/-/_/g');do
    find $closed_dir -name "pFind-Filtered.spectra"|xargs cat|grep $sample > $closed_dir/$sample.spectra.tmp
    cat $closed_dir/closed.res.header.txt $closed_dir/$sample.spectra.tmp > $closed_dir/$sample.spectra
    find $open_dir -name "pFind-Filtered.spectra"|xargs cat|grep $sample > $open_dir/$sample.spectra.tmp
    cat $open_dir/open.res.header.txt $open_dir/$sample.spectra.tmp > $open_dir/$sample.spectra
    rm -rf $closed_dir/$sample.spectra.tmp $open_dir/$sample.spectra.tmp
done