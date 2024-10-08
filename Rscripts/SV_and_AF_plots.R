#vvvvvvvvv with tsv
title <- "adotto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed_oc_ready.vcf.tsv"
tsv_file <- read.csv(title,sep="\t",skip=5)

# Use All Mappings column for gene name; between first and second :
tsv_Gene <- strsplit(tsv_file$Gene)
AF <- list(tsv_file$AC/tsv_file$AN)







#vvvvvvv with excel file

vcf_file <- read_xlsx("addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed_oc_ready.vcf.xlsx", "Gene", skip=1)
variant_file <- read_xlsx("addoto_variants.grch38.sqoff.chr1.above50bp.truvari_collapsed_oc_ready.vcf.xlsx", "Variant", skip=1)

coding <- data.frame(Genes=vcf_file$Gene,Number_of_Variants=vcf_file$`Number of Coding Variants`)
noncoding <- data.frame(name=vcf_file$Gene,value=vcf_file$`Number of Noncoding Variants`)

result <- variant_file %>% 
  filter(Gene>0)

x <- result$Gene
y <- result$AC/result$AN

#install.packages("ggplot2")
library(ggplot2)
library(dplyr)

test_data <- aes(x=x,y=y)
test_graph <- ggplot(x, test_data)+geom_violin(aes(colour="#1268FF"),alpha=0.3)
print(test_graph)

new_coding <- aes(x=reorder(Genes,-Number_of_Variants),y=Number_of_Variants)

coding_graph <- ggplot(head(coding,10), new_coding)+geom_bar(stat='identity')+
  xlab("") + ylab("Number of Variants")+ ggtitle("Top 10 Coding Variants in Genes for Chr1")

print(coding_graph)
