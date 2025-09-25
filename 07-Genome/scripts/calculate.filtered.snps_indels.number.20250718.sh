source activate biotools

INPUT_VCF="/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter/human_brain_21pcw_filtered_pass.vcf.gz"
OUTPUT_DIR="/home/user/data3/lit/project/sORFs/07-Genome/processed/vcf_filter"
SAMPLE_NAME="human_brain_21pcw"

gatk SelectVariants \
    -V ${INPUT_VCF} \
    -select-type SNP \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_snps_filtered_pass.vcf.gz

gatk SelectVariants \
    -V ${INPUT_VCF} \
    -select-type INDEL \
    -O ${OUTPUT_DIR}/${SAMPLE_NAME}_indels_filtered_pass.vcf.gz

echo "过滤后的SNP数量:" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report_add.txt
bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_snps_filtered_pass.vcf.gz | wc -l >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report_add.txt
echo "过滤后的INDEL数量:" >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report_add.txt
bcftools view -H ${OUTPUT_DIR}/${SAMPLE_NAME}_indels_filtered_pass.vcf.gz | wc -l >> ${OUTPUT_DIR}/${SAMPLE_NAME}_filtering_report_add.txt