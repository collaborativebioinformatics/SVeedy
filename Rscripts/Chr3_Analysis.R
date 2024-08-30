

# Read tsv file
Chr03 <- read.table(file='adotto_variants.grch38.sqoff.chr3.above50bp.truvari_collapsed_oc_ready.vcf.tsv', skip = 6, sep = '\t', header =  FALSE, fill = TRUE)
header <- read.table(file='adotto_variants.grch38.sqoff.chr3.above50bp.truvari_collapsed_oc_ready.vcf.tsv', sep = '\t', skip = 2, nrows =1,header =  TRUE, fill = TRUE)
colnames(Chr03) <- unlist(header)

# Count clinical significant variants
nrow(Chr03[Chr03$`Clinical Significance` == "Pathogenic",]) # returns 3
nrow(Chr03[Chr03$`Clinical Significance` == "Likely pathogenic",]) # returns 0

# Subset Pathogenic variants
Chr3_Pathogenic <- Chr03[Chr03$`Clinical Significance` == 'Pathogenic', ]

# Write data frame into csv file
write.csv(Chr3_Pathogenic, "Chr3_Pathogenic.csv")
