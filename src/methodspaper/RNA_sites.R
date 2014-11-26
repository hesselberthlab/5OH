# Calculate rRNA, mRNA, tRNA, etc distributions for all samples from GFF intersection

library(ggplot2)
library(reshape2)
library(dplyr)

# Return table with sample, site, count
get_RNA_site_counts <- function(sample.name, fileprefix) {
  intersectfile <- paste("rgraphs/", fileprefix, ".gffintersect.tab", sep="")
  all_data <- read.table(intersectfile, header=TRUE)
  data <- all_data %>%
    select(site, count)
  site.counts <- data[!duplicated(data), ] %>%
    mutate(sample = sample.name) %>%
    select(sample, site, count)
  return (site.counts)
}

# Directory containing results
setwd("~/projects/5OH/results/methodspaper/")

# Useful library prefixes
sample_libraries = c("SP5", "SP6", "SP7", "SP8", "SP9", "SP10", "SP11", "SP12", "SP13", "SP14", "SP15", "SP16", "SP17", "SP18", "SP19", "SP27")
noSAP_libraries = c("SP5", "SP6", "SP7", "SP8", "SP13", "SP14", "SP15", "SP27")
SAP_libraries = c("SP9", "SP10", "SP11", "SP12", "SP16", "SP17", "SP18", "SP19")
control_libraries = c("SP20", "SP28", "SP29")
all_libraries = c("SP5", "SP6", "SP7", "SP8", "SP9", "SP10", "SP11", "SP12", "SP13", "SP14", "SP15", "SP16", "SP17", "SP18", "SP19", "SP27", "SP20", "SP28", "SP29")

# Alignment specifics
assembly <- "assembly-sacCer1"
alignment <- "align-all"

# Create metadata table of all samples and their RNA site counts
sitemetadata <- data.frame()
for (sample.name in all_libraries) {
  fileprefix <- paste(sample.name, assembly, alignment, sep=".")
  data <- get_RNA_site_counts(sample.name, fileprefix)
  sitemetadata <- rbind(sitemetadata, data)
}

rRNA5.8 <- "chr12:455571"
subset(sitemetadata, site == rRNA5.8)

SNR18 <- "chr1:142368"
subset(sitemetadata, site == SNR18)