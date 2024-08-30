#!/bin/Rscript
options(stringsAsFactors = F)

args = commandArgs(TRUE)

allCallsfl=args[1] # SURVIVOR Stats for file with all calls
above50fl=args[2] # SURVIVOR Stats for file with calls with < 50bp length filtered out
truvarifl=args[3] # SURVIVOR Stats for file with Truvari collapsed/merged SV calls 
workdir=args[4] # Path to input and output files

library("tidyverse")
library("ggplot2")
library("scales") 

setwd(workdir)

all <- read.table(allCallsfl, header = T, sep = "\t")
abv50bp <- read.table(above50fl, header = T, sep = "\t")
abv50bpTruvari <- read.table(truvarifl, header = T, sep = "\t")

all$condition <- "All"
abv50bp$condition <- ">50bp"
abv50bpTruvari$condition <- ">50bp and Merged"

df <- rbind.data.frame(all, abv50bp, abv50bpTruvari)
df2 <- df[,c(8, 1:6)]

dfFinal <- df2 %>% pivot_longer(
  cols = c("Del", "Dup", "Inv", "INS", "TRA"),
  names_to = "SVtype",
  values_to = "Count"
)

conditionF <- factor(dfFinal$condition, levels = c('All', '>50bp', '>50bp and Merged'))

png(filename = "SURVIVOR-stats.png", res=300, height=4.5, width=6.5, units="in")
ggplot(dfFinal) +
  geom_bar(aes(x = conditionF, y = Count, fill = SVtype),
           position = "dodge",
           stat = "identity") +
  facet_grid(~ factor(dfFinal$Len, levels=c('0-50bp','50-100bp','100-1000bp','1000-10000bp','10000+bp'))) + theme_bw() +
  theme(axis.text.x = element_text(angle=60,hjust = 1), axis.title.x = element_blank()) +
  scale_y_continuous(labels = scales::comma)

dev.off()
