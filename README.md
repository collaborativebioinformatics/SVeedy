# SVeedy
# Contributors: 
1. Jon Moller (CARD, NIA, NIH, Bethesda, Maryland)
2. Sam Stroupe (Texas A&M University)
3. Sarah Fross (Texas A&M University)
4. Shaghayegh Beheshti (Baylor College of Medicine)
5. Pankhuri Wanjari (University of Chicago- Department of Pathology)
6. Nha Van Huynh (University of Alabama at Birmingham)
![image](https://github.com/user-attachments/assets/5bf7b126-30c8-48cd-8d2f-f0049b398b5e)
# Structural Variant Analysis Pipeline

## Overview
Structural variants (SVs) represent deviations from a reference genome sequence, typically spanning more than 50 base pairs (bps). These variations can have significant implications for understanding genetic diversity and the mechanisms underlying various phenotypes. This project aims to develop a robust pipeline for detecting and cataloging identical SVs across different samples and databases, ultimately linking them to specific phenotypes.
# Structural Variants in Human Genomes
![image](https://github.com/user-attachments/assets/e8a32170-29e0-4fdd-8700-28db0f4be256)
Larger structural variants are present among human genomes. For example, human chromosomes can have:
- **Missing segments** (*deletion variants*)
- **Duplicated segments** (*duplication variants*)
- **Inverted segments** (*inversion variants*)
- **Added segments** (*insertion variants*)
- **Segments transferred from other chromosomes** (*translocation variants*)

*Source: [NIH Human Genomic Fact Sheet](https://www.genome.gov/about-genomics/fact-sheets/Genomic-Research) (Last updated: February 1, 2023)*
## Objectives
The primary goal of this study is to identify and analyze SVs in novel and known genes, as well as established population SVs, to uncover new biological processes and associations. By cross-referencing SVs with phenotypic data, this pipeline seeks to establish a more comprehensive understanding of genotype-phenotype correlations.

## Methodology
- **SV Detection:** The pipeline will accurately identify SVs across multiple datasets, ensuring consistency and reliability in detecting known and novel variants.
- **Phenotype Association:** Each identified SV will be linked to phenotypic data, allowing for the correlation of specific genetic variations with particular traits or diseases.
- **VCF File Output:** The results will be condensed into a variant calling format (VCF) file, summarizing the detected SVs and their associated phenotypes. Users can then input a patient ID to retrieve potential phenotypic outcomes based on the identified SVs.

## Case Study
As an example, a previous study identified 11 SV loci associated with an increased risk for obesity, with an Odds Ratio exceeding 25% ([DOI: 10.1371/journal.pone.0058048](https://doi.org/10.1371/journal.pone.0058048)). This project aims to build upon such findings by extending the analysis to a broader set of SVs and phenotypes, facilitating the discovery of novel genetic contributors to complex traits.

# Workflow
![SVeedy - Wed](https://github.com/user-attachments/assets/54dfa671-0fee-4056-91d5-cfdfa21c287e)

![SVeedy - Fri](https://github.com/user-attachments/assets/68e8a3c2-73af-4aa8-9755-7bbc30abb8ee)


# Validation

We are validating our pipeline on the Project Adotto assembly-based variant calls from the GIAB tandem repeat benchmark (https://zenodo.org/records/6975244), beginning with SV calls in chromosome 1 (either insertions or deletions). Upon SV filtering steps, we went from 194,098 SVs to 55,905 SVs (remove those under 50 bp in length) and then 29,026 SVs (truvari collapse function keeping most common allele in each cluster).

![image](https://github.com/user-attachments/assets/5ec8c25c-7c7d-412d-a323-93d087d679e3)
Above: Raw SV counts after filtering (>50 bp SVs) and merging (truvari collapse) steps of pipeline for Project Adotto chromosome 1 variants.

# Demographic features of the study population 
<img width="500" alt="image" src="https://github.com/user-attachments/assets/206c3939-ada4-445c-ac75-cf93471ebab9">
<img width="467" alt="image" src="https://github.com/user-attachments/assets/4f470f03-b84c-49cc-bba4-001bc855eff8">

# Gene Analysis of Chromosome 1

In our analysis of chromosome 1 using the Adotto dataset, we identified genes with the most prevalent allele frequencies across different populations. These allele frequencies, including those for structural variants, were sourced from the gnomAD dataset. This analysis highlights genes that show significant variation in allele frequencies among American, Ashkenazi Jewish, East Asian, Finnish, Non-Fin European, and Other populations.

For example, the gene **NFASC**, which is involved in neurodevelopmental disorders with central and peripheral motor dysfunction (MIM 609145), shows notable structural variants in the East Asian ancestry. The prevalence of structural variants of NFASC in this population underscores the importance of understanding population-specific genetic variations, which can inform physicians and researchers about potential genetic risk factors and guide future studies. Our tool with continued research into these population-specific variants is essential for advancing personalized medicine and improving genetic counseling.
![image](https://github.com/user-attachments/assets/fc0d8dc3-c494-4992-af07-9e5c413db348)

# Conclusion and Future Directions

