library(dplyr)      ## v1.1.4
library(tidyr)      ## v1.3.1
library(argparse)   ## v2.2.3


parser <- ArgumentParser()
parser$add_argument("-t", "--tsv", help="OpenCRAVAT report tsv output")
args <- parser$parse_args()



##### Read in the OpenCRAVAT tsv report and skip the header lines #####
report <- read.delim(args$tsv, skip = 5)
#report <- read.delim("CravatReportTest.tsv", skip = 5)

# Keep only the variants that have a known clinical significance
disease_associated <- report %>%
  filter(Clinical.Significance %in% c('pathogenic', 'likely pathogenic')) %>%
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
  
  # Set margins
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
      "Sample ", sample_id, " likely has ", variant$Disease.Names,
      " due to the identified structural variant at the genomic coordinate ", variant$Chrom,
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

  
