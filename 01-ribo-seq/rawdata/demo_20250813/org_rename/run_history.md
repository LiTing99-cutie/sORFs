for f in ../MJ20250804226-ZX-D-250801018-李春琼-纯文库-18个样本/rawdata/*R1.raw.fastq.gz; do
    newname=$(basename "$f" | sed 's/\.R1\.raw\.fastq\.gz$/.fq.gz/')
    ln -s "$(realpath "$f")" "./$newname"
done