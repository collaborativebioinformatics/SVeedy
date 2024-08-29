#!/bin/bash
# sveedy preprocessing bash script optimized for Biowulf
# run in interactive node with 32GB RAM and 32 cores
sinteractive --mem=32 --cpus-per-task=32
# following dependencies: samtools, bcftools, truvari, SURVIVOR
# first convert bcf to vcf
# change folders as necessary
BASE_DIR=/data/mollerabg/BCM_HGSC_hackathon_2024/
cd ${BASE_DIR}/Jon_chr1_tests_082924
# use bcftools to convert bcf to vcf
# bcftools/samtools 1.19
module load bcftools
bcftools convert -O v -o addoto_variants.grch38.sqoff.chr1.vcf ../adotto_variants.grch38.sqoff.chr1.bcf.gz
# filter out SVs under 50 bp
bcftools view -i "SVLEN > 50" addoto_variants.grch38.sqoff.chr1.vcf -o addoto_variants.grch38.sqoff.chr1.above50bp.vcf
# collapse variants with truvari
# Truvari 4.3.0
# make sure truvari is installed with "pip install truvari"
# index and bgzip variants first
bgzip -c addoto_variants.grch38.sqoff.chr1.above50bp.vcf > addoto_variants.grch38.sqoff.chr1.above50bp.vcf.gz
tabix -p vcf addoto_variants.grch38.sqoff.chr1.above50bp.vcf.gz
# then run truvari
module load python
truvari collapse -k common -i addoto_variants.grch38.sqoff.chr1.above50bp.vcf.gz -o addoto_variants.grch38.sqoff.chr1.above50bp.truvari_raw.vcf -c addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed.vcf
# get statistics of SVs at each stage to share
# use SURVIVOR 1.0.7 (from Sedlazcek lab)
module load survivor
# initial VCF
SURVIVOR stats addoto_variants.grch38.sqoff.chr1.vcf -1 -1 -1 addoto_variants.grch38.sqoff.chr1.survivor_stats
# after filtering SVs under 50 bp
SURVIVOR stats addoto_variants.grch38.sqoff.chr1.above50bp.vcf -1 -1 -1 addoto_variants.grch38.sqoff.chr1.above50bp.survivor_stats
# after truvari common collapse
SURVIVOR stats addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed.vcf -1 -1 -1 addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed.survivor_stats
# now run open-cravat annotation on collapsed SVs
module load opencravat
# remove problematic header line
bcftools annotate -x "FILTER/COV" addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed.vcf > addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed_oc_ready.vcf
oc run addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed_oc_ready.vcf -l hg38 -a gnomad clinvar -t text excel -d addoto_variants_chr1_above50bp_truvari_collapsed_oc
