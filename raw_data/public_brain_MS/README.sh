source activate fragpipe
cd /rd1/user/lit/project/sORFs/raw_data/public_brain_MS
[ -f HLA/PXD019643/ms2_count.txt ] && rm -rf HLA/PXD019643/ms2_count.txt
for raw in $(find HLA/PXD019643/ -name "*.raw");do
ms2_count=$(mono /rd1/user/lit/software/ThermoRawFileParser.exe -i $raw -s | grep -c "<spectrum ")
echo -e "$raw\t$ms2_count"
done >> HLA/PXD019643/ms2_count.txt

find HLA/PXD019643/ -name "*.raw" | parallel -j 8 'ms2_count=$(mono /rd1/user/lit/software/ThermoRawFileParser.exe -i {} -s | grep -c "<spectrum "); echo -e "{}\t$ms2_count"' >> HLA/PXD019643/ms2_count.txt

#32542
mono /rd1/user/lit/software/ThermoRawFileParser.exe -i non-HLA/PXD035950/NGN2_87922.raw -s | grep -c "<spectrum "
#32542
grep -c "<spectrum " non-HLA/PXD035950/NGN2_87925_calibrated.mzML
less non-HLA/PXD035950/NGN2_87925_calibrated.mzML | grep 'ms level' | grep -c 'value="2"'
#32563
grep -c "<spectrum " non-HLA/PXD035950/NGN2_87925_uncalibrated.mzML
#32563
less non-HLA/PXD035950/NGN2_87925_uncalibrated.mzML | grep 'ms level' | grep -c 'value="2"'
#36302
mono /rd1/user/lit/software/ThermoRawFileParser.exe -i non-HLA/PXD035950/NGN2_87922.raw -s | grep -c "<spectrum "
# value="2" name="ms level"
#33664
mono /rd1/user/lit/software/ThermoRawFileParser.exe -i non-HLA/PXD035950/NGN2_87922.raw -s | grep 'ms level' | grep -c 'value="2"'
mono /rd1/user/lit/software/ThermoRawFileParser.exe -i non-HLA/PXD035950/NGN2_87922.raw -s > non-HLA/PXD035950/tmp/NGN2_87922.mzML
grep 33664 non-HLA/PXD035950/tmp/NGN2_87922.mzML
grep 36302 non-HLA/PXD035950/tmp/NGN2_87922.mzML
#33664
#在windows上进行msconvert，并且查找level" value="2"的行数

# 还没运行完，但是先不运行
[ -f HLA/PXD019643/ms2_count.txt ] && rm -rf HLA/PXD019643/ms2_count.txt
for raw in $(find HLA/PXD019643/ -name "*.raw");do
ms2_count=$(mono /rd1/user/lit/software/ThermoRawFileParser.exe -i $raw -s | grep 'ms level' | grep -c 'value="2"')
echo -e "$raw\t$ms2_count"
done >> HLA/PXD019643/ms2_count.txt