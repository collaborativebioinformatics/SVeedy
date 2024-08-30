# SVeedy
![image](https://github.com/user-attachments/assets/5bf7b126-30c8-48cd-8d2f-f0049b398b5e)

# Structural Variant Analysis Pipeline

## Overview & Objective
Structural variants (SVs) represent deviations from a reference genome sequence, typically spanning more than 50 base pairs (bps). These variations can have significant implications for understanding genetic diversity and the mechanisms underlying various phenotypes. This project aims to develop a robust pipeline for detecting and cataloging identical SVs across different samples and databases, ultimately linking them to specific phenotypes. The primary goal of this study is to identify and analyze SVs in novel and known genes, as well as established population SVs, to uncover new biological processes and associations. By cross-referencing SVs with phenotypic data, this pipeline seeks to establish a more comprehensive understanding of genotype-phenotype correlations.

# Structural Variants in Human Genomes
![image](https://github.com/user-attachments/assets/e8a32170-29e0-4fdd-8700-28db0f4be256)
Larger structural variants are present among human genomes. For example, human chromosomes can have:
- **Missing segments** (*deletion variants*)
- **Duplicated segments** (*duplication variants*)
- **Inverted segments** (*inversion variants*)
- **Added segments** (*insertion variants*)
- **Segments transferred from other chromosomes** (*translocation variants*)

*Source: [NIH Human Genomic Fact Sheet](https://www.genome.gov/about-genomics/educational-resources/fact-sheets/human-genomic-variation) (Last updated: February 1, 2023)*

## Methodology
- **SV Detection:** The pipeline will accurately identify SVs across multiple datasets, ensuring consistency and reliability in detecting known and novel variants.
- **Phenotype Association:** Each identified SV will be linked to phenotypic data, allowing for the correlation of specific genetic variations with particular traits or diseases.
- **VCF File Output:** The results will be condensed into a variant calling format (VCF) file, summarizing the detected SVs and their associated phenotypes. Users can then input a patient ID to retrieve potential phenotypic outcomes based on the identified SVs.

We gained access to a collection of VCFs created to find Tandem Repeats (TRs) (English et. al 2024)  from a collection of 86 haplotypes accumulated from Garg et. al 2020, Ebert et. al 2021, Jarvis et. al 2022, and Wang et. al 2022.

Then, we converted these BCF files to VCF files and filtered out SVs under 50 bp:

`bcftools convert -O v -o file.vcf file.bcf.gz`

`bcftools view -i "SVLEN > 50" file.vcf -o above50bp.vcf`

One thing we need to consider is the cutoff percentage of similarity of the SVs. For example, if obesity has been identified to have a 100bp SV compared to the population reference, we would want to determine if 80 out of the 100 bps (80%) are the same. This is one of the goals of the Truvari software (English et. al 2022) and they explain that although a similar SV may be present in different samples, it can be in different loci along the genome. They also explain that over-filtering using the collapse command may remove regions of the SV. We also need to be wary of defining the range of each SV and whether we will add, for example, a buffer of 5bp on either side of a SV to take unique alleles into account.

By using the Truvari software, we were able to collapse the SVs:

`truvari collapse -k common -i above50bp.vcf.gz -o truvari_raw.vcf -c truvari_collapsed.vcf`

At this point, we could run SURVIVOR 1.0.7 [DOI: doi:10.12688/f1000research.12516.1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5668921/) for analyzing the VCF data:

`SURVIVOR stats file.vcf -1 -1 -1 survivor_stats`

Finally, we used OpenCRAVAT (Pagal et. al 2019) to annotate the VCF file using the ClinVar and gnomAD databases and hg38 genome reference. We removed a problematic line of the VCF header (FILTER/COV) before running the input collapsed SV VCF through open-cravat.

`oc run above50bp.truvari_collapsed_oc_ready.vcf -l hg38 -a gnomad clinvar -t text excel -d opencravat_output --debug`

We also added in cosmic, gnomad_gene, clinvar_acmg annotators

`oc run above50bp.truvari_collapsed_oc_ready.vcf -l hg38 -a gnomad gnomad_gene clinvar clinvar_acmg cosmic cosmic_gene -t text excel -d opencravat_output2 --debug`

## Case Study
As an example, a previous study identified 11 SV loci associated with an increased risk for obesity, with an Odds Ratio exceeding 25% ([DOI: 10.1371/journal.pone.0058048](https://doi.org/10.1371/journal.pone.0058048)). This project aims to build upon such findings by extending the analysis to a broader set of SVs and phenotypes, facilitating the discovery of novel genetic contributors to complex traits.


# Workflow

![pankhuri's flowchart](https://github.com/user-attachments/assets/cc6ce172-667b-4663-9873-e2f257e02470)

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

![10_coding_chr1](https://github.com/user-attachments/assets/0448ad38-1404-44f6-9e0b-26bc37f2dd37)


# Annotating SVs with ClinVar Annotations

We have successfully validated our pipeline by gathering all structural variants (SVs) from the Adotto database and combining them with ClinVar data. This analysis led to the identification of **three structural variants** classified as pathogenic in ClinVar, affecting a total of **eight individuals** across the dataset.


The figure below illustrates the results of our analysis, including a table that provides detailed information about the identified pathogenic structural variants:
<img width="900" alt="All_ClinVar_hits" src="https://github.com/user-attachments/assets/dc256595-5b74-4b0b-8a07-acd969307d4c">
- **ClinVar ID**: The unique identifier for each variant in the ClinVar database.
- **Disease Name**: The associated disease or condition, along with corresponding Human Phenotype Ontology (HPO) terms.
- **Clinical Significance**: The clinical classification of the variant, which in this analysis, is labeled as pathogenic.
![image](https://github.com/user-attachments/assets/d6fa1ce4-e1d1-4753-b615-6014300cc735)

This table serves as a reference for the identified structural variants and their clinical implications, based on the combined analysis of the Adotto database and ClinVar data.

## User-Friendly Tool for PDF Generation

To facilitate the use of this information, we have developed an additional tool that converts the data into a user-friendly PDF output. This PDF includes the sample name of each individual and the predicted diagnosis based on the known ClinVar phenotypes. It provides a comprehensive report detailing each variant, including all relevant information.

An example of this PDF report is shown below:

![image](https://github.com/user-attachments/assets/e7d56956-101c-4171-97cc-5bf589d7fcaa)

This tool streamlines the interpretation of SV data and enhances the ease of accessing detailed diagnostic information, making it an invaluable resource for clinical research and patient care.

# Future Directions and Conclusion

**Future Directions:**

- **Enhanced Integration**: Expand the pipeline by incorporating additional databases and resources to improve the breadth and depth of structural variant (SV) annotations.

- **Binder Implementation**: Integrate the entire workflow into a Binder environment for increased accessibility and reproducibility of analyses.

- **Streamlining and Optimization**: Focus on streamlining the pipeline to ensure faster and more efficient processing of SV data, including optimizing data flow and reducing computational overhead.

- **Methodological Expansion**: Explore and implement additional methodologies within the workflow to enhance the accuracy and scope of SV detection and annotation.

**Conclusion:**

The development of this pipeline represents a significant advancement in the annotation of structural variants (SVs). By combining gnomAD allele frequencies and ClinVar clinical data, our tool facilitates a more straightforward and efficient approach to detecting and analyzing SVs in patient sequences. The integration of phenotypic information with clinical and larger dataset sources enhances the tool's utility in patient care, leading to more informed predictions and better clinical decision-making. The streamlined design and future expansions aim to set a new standard for bioinformatics workflows in precision medicine.

# Contributors: 
1. Jon Moller (CARD, NIA, NIH, Bethesda, Maryland)
2. Sam Stroupe (Texas A&M University)
3. Sarah Fross (Texas A&M University)
4. Shaghayegh Beheshti (Baylor College of Medicine)
5. Pankhuri Wanjari (University of Chicago- Department of Pathology)
6. Nha Van Huynh (University of Alabama at Birmingham)
