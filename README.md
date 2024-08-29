# SVeedy
# Contributors: 
1. Jon Moller (CARD, NIA, NIH, Bethesda, Maryland)
2. Sam Stroupe (Texas A&M University)
3. Sarah Fross (Texas A&M University)
4. Shaghayegh Beheshti (Baylor College of Medicine)
5. Pankhuri Wanjari (University of Chicago- Department of Pathology)
6. Nha Van Huynh (University of Alabama at Birmingham)
7. Ali Saadat

![image](https://github.com/user-attachments/assets/5bf7b126-30c8-48cd-8d2f-f0049b398b5e)

# Structural Variant Analysis Pipeline

## Overview
Structural variants (SVs) represent deviations from a reference genome sequence, typically spanning more than 50 base pairs (bps). These variations can have significant implications for understanding genetic diversity and the mechanisms underlying various phenotypes. This project aims to develop a robust pipeline for detecting and cataloging identical SVs across different samples and databases, ultimately linking them to specific phenotypes.

## Objectives
The primary goal of this study is to identify and analyze SVs in novel and known genes, as well as established population SVs, to uncover new biological processes and associations. By cross-referencing SVs with phenotypic data, this pipeline seeks to establish a more comprehensive understanding of genotype-phenotype correlations.

## Methodology
- **SV Detection:** The pipeline will accurately identify SVs across multiple datasets, ensuring consistency and reliability in detecting known and novel variants.
- **Phenotype Association:** Each identified SV will be linked to phenotypic data, allowing for the correlation of specific genetic variations with particular traits or diseases.
- **VCF File Output:** The results will be condensed into a variant calling format (VCF) file, summarizing the detected SVs and their associated phenotypes. Users can then input a patient ID to retrieve potential phenotypic outcomes based on the identified SVs.

## Case Study
As an example, a previous study identified 11 SV loci associated with an increased risk for obesity, with an Odds Ratio exceeding 25% ([DOI: 10.1371/journal.pone.0058048](https://doi.org/10.1371/journal.pone.0058048)). This project aims to build upon such findings by extending the analysis to a broader set of SVs and phenotypes, facilitating the discovery of novel genetic contributors to complex traits.

# Workflow
![image](https://github.com/user-attachments/assets/ea484571-3640-42b6-acdc-3926434a07ae)

