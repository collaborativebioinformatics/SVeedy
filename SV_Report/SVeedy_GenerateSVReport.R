library(dplyr)      ## v1.1.4
library(tidyr)      ## v1.3.1
library(argparse)   ## v2.2.3


parser <- ArgumentParser()
parser$add_argument("-t", "--tsv", action = "append",
                    help="OpenCRAVAT report tsv file", required = TRUE)
args <- parser$parse_args()


##### Read in the OpenCRAVAT tsv report and skip the header lines #####
# store it as data frames in a list
#report <- read.delim(args$tsv, skip = 5)
# report1 <- read.delim("adotto_chr1.tsv", skip = 5)
# report2 <- read.delim("adotto_chr2.tsv", skip = 5)
# report3 <- read.delim("adotto_chr3.tsv", skip = 5)
# list <- c("adotto_chr1.tsv","adotto_chr2.tsv","adotto_chr3.tsv")

# Initialize an empty list to store all of the tsv files
data <- list()
for (file in args$tsv) {
  report <- read.delim(file, skip = 5)
  data[[length(data) + 1]] <- report
}
# for (file in list) {
#   report <- read.delim(file, skip = 5)
#   data[[length(data) + 1]] <- report
# }

# Merge all data frames into one report
report <- do.call(rbind, data)

###### Analyze the data ############
# Make plots on general trends of the data
## GET CODE FROM GROUP ##


##### Identify Variants of Clinical interest in the data base ##

# Keep only the variants that have a known clinical significance
disease_associated <- report %>%
  filter(Clinical.Significance %in% c('Pathogenic', 'Likely Pathogenic', 'VUS')) %>%
  separate_rows(Samples, sep = ";")

# Group by Sample and summarize the associated variants
variants_per_sample <- disease_associated %>%
  group_by(Samples) %>%
  summarize(Associated_Variants = list(data.frame(Chrom, Position, Clinical.Significance,
                                                  Disease.Names, Disease.Ref.Nums)))

#### Generate a report based on the file ####

# Open PDF device for the multi-page report
pdf(file = "Genetic Report of Diseases caused by SVs.pdf")

# Loop through each sample to add a new page to the PDF report
for (i in 1:nrow(variants_per_sample)) {
  sample_id <- variants_per_sample$Samples[i]
  associated_variants <- variants_per_sample$Associated_Variants[[i]]
  
  # Start a new page for each sample
  plot.new()
  
  # Write a title to the new page
  title(main = paste("Genetic Report for Sample", sample_id))
  
  # Set margins for graphs if needed
  #par(mar = c(4, 4, 4, 4))
  
  # Set initial y position for text
  y_position <- 1
  
  # Set line height and max width for text wrapping
  line_height <- 0.05
  max_width <- 80  # Adjust this width as necessary
  
  # Write the paragraph for each variant
  for (j in 1:nrow(associated_variants)) {
    variant <- associated_variants[j, ]
   
  # Example paragraph that can be written into a generated report that
    # inserts information taken from the OpenCRAVAT report 
     paragraph <- strwrap(paste(
      "Sample ", sample_id, " is at risk for ", variant$Disease.Names,
      " due to an identified structural variant at the genomic coordinate ", variant$Chrom,
      ":", variant$Position, ". The identified structural variant 
      in the genome of this individual is ",
      variant$Clinical.Significance, " based on previous research. 
      More information about this genetic variant can be found in the 
      ClinVar database using the ClinVar reference ID ",
      variant$Disease.Ref.Nums, ".", sep = ""
    ), width = max_width)
    
    # Adjust y position for each line of wrapped text
    for (line in paragraph) {
      if (y_position - line_height < 0) {
        plot.new()  # Start a new page if there's no more space on the current page
        y_position <- 1
      }
      # 
      text(x = 0, y = y_position, labels = line, pos = 4, cex = 0.8)
      y_position <- y_position - line_height
    }
  }
}

# Close the PDF device
dev.off()

  
