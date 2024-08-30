#!/usr/bin/env bash

#SBATCH --cpus-per-task=32
#SBATCH --mem=80g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=8:00:00
#SBATCH --gres=lscratch:20

# sveedy preprocessing bash script optimized for Biowulf
# slurm parameters above
# run list of inputs as array job
# run in interactive node with 32GB RAM and 32 cores
# comment out sinteractive --mem=32 --cpus-per-task=32
# following dependencies: samtools, bcftools, truvari, SURVIVOR
# first convert bcf to vcf
# change folders as necessary
BASE_DIR=/data/mollerabg/BCM_HGSC_hackathon_2024/
cd ${BASE_DIR}
# eventually create command line parameter for FILELIST
FILELIST=${BASE_DIR}/project_addoto_bcfs_083024.txt
# set SAMPLE_ID for particular job in array based on task id
N=${SLURM_ARRAY_TASK_ID}
SAMPLE_ID=$(sed -n ${N}p ${FILELIST} | sed 's/.bcf.gz//g')
# convert periods to underscores
FOLDER_NAME=$(echo ${SAMPLE_ID} | sed 's/\./_/g')
# important info for the slurm logs
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FOLDER_NAME ${FOLDER_NAME}"

# make separate directories for each set of outputs
mkdir -p ${BASE_DIR}/${FOLDER_NAME}

# use bcftools to convert bcf to vcf
# bcftools/samtools 1.19
module load bcftools
bcftools convert -O v -o ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.vcf ${BASE_DIR}/${SAMPLE_ID}.bcf.gz
# filter out SVs under 50 bp
bcftools view -i "SVLEN > 50" ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.vcf -o ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.vcf
# collapse variants with truvari
# Truvari 4.3.0
# make sure truvari is installed with "pip install truvari"
# index and bgzip variants first
bgzip -c ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.vcf > ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.vcf.gz
tabix -p vcf ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.vcf.gz
# then run truvari
module load python
truvari collapse -k common -i ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.vcf.gz -o ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_raw.vcf -c ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed.vcf
# get statistics of SVs at each stage to share
# use SURVIVOR 1.0.7 (from Sedlazcek lab)
module load survivor
# initial VCF
SURVIVOR stats ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.vcf -1 -1 -1 ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.survivor_stats
# after filtering SVs under 50 bp
SURVIVOR stats ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.vcf -1 -1 -1 ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.survivor_stats
# after truvari common collapse
SURVIVOR stats ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed.vcf -1 -1 -1 ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed.survivor_stats
# now run open-cravat annotation on collapsed SVs
# use open-cravat 2.4.2
module load opencravat
# remove problematic header line
bcftools annotate -x "FILTER/COV" ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed.vcf > ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed_oc_ready.vcf
oc run ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed_oc_ready.vcf -l hg38 -a gnomad clinvar -t text excel -d ${BASE_DIR}/${FOLDER_NAME}/opencravat_output --debug
# add in cosmic, gnomad_gene, clinvar_acmg annotators
oc run ${BASE_DIR}/${FOLDER_NAME}/${SAMPLE_ID}.above50bp.truvari_collapsed_oc_ready.vcf -l hg38 -a gnomad gnomad_gene clinvar clinvar_acmg cosmic cosmic_gene -t text excel -d ${BASE_DIR}/${FOLDER_NAME}/opencravat_output2 --debug



