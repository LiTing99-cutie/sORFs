set -eo pipefail

source activate base

output_path=$1
mkdir -p $output_path && cd $output_path
Ribocode_bam=$2
PRICE_bam=$3
RibORF_sam=$4

mkdir -p $output_path/RiboCode && cd $output_path/RiboCode
mkdir -p $output_path/PRICE && cd $output_path/PRICE