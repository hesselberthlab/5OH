# Calculate rRNA, mRNA, tRNA, etc distributions for all samples from GFF intersection

library(ggplot2)
library(reshape2)
library(dplyr)

# Calculate total number of reads that went into GFF intersection
calculate_reads <- function(fileprefix) {
  pos_bg <- paste(fileprefix, ".strand.pos.CPMs.bg", sep="")
  neg_bg <- paste(fileprefix, ".strand.neg.CPMs.bg", sep="")
  
  pos_counts <- read.table(pos_bg, col.names = c("chrom", "start", "stop", "count")) %>%
    summarize(total = sum(count))
  
  neg_counts <- read.table(neg_bg, col.names = c("chrom", "start", "stop", "count")) %>%
    summarize(total = sum(count))
  
  total <- pos_counts$total + neg_counts$total
  
  return (total)
}

# Return a data table with total counts of RNA subspecies, pct of total counts, and sample name
get_RNA_species_count <- function(sample.name, fileprefix, RNAtypes) {
  total_reads <- calculate_reads(fileprefix)
  intersectfile <- paste("rgraphs/", fileprefix, ".gffintersect.tab", sep="")
  all_data <- read.table(intersectfile, header=TRUE)
  data <- all_data %>%
    select(chr, start, count, cat)
  species.counts <- subset(data[!duplicated(data), ], cat %in% RNAtypes) %>%
    group_by(cat) %>%
    summarize(counts = sum(count)) 
  data <- rbind(species.counts, data.frame(cat="total", counts=total_reads)) %>%
    mutate(pct = counts/total_reads, sample = sample.name)
  return (data)
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

# Create metadata table of all samples and their RNA species counts
RNAtypes <- c("rRNA", "mRNA", "5UTR", "3UTR", "tRNA")
metadata <- data.frame()
for (sample.name in all_libraries) {
  fileprefix <- paste(sample.name, assembly, alignment, sep=".")
  data <- get_RNA_species_count(sample.name, fileprefix, RNAtypes)
  metadata <- rbind(metadata, data)
}


# Summary plots of alignment data
ggplot(subset(metadata, cat != "total"), aes(x=cat, y=pct, fill=sample)) +
  geom_bar(stat="identity", position="dodge")

# Calculating percentages for subsets of data
pcttable <- data.frame(subset(metadata, sample %in% noSAP_libraries)) %>%
  select(sample, cat, counts) %>%
  group_by(cat) %>%
  summarize(total.counts = sum(counts))

total <- subset(pcttable, cat == "total")$total.counts

pcttable <- pcttable %>%
  mutate(pct = total.counts/total)


# Investigation of RiboZero sample's high tRNA content
for (sample.name in c("SP28")) {
  fileprefix <- paste(sample.name, assembly, alignment, sep=".")
  intersectfile <- paste("rgraphs/", fileprefix, ".gffintersect.tab", sep="")
  all_data <- read.table(intersectfile, header=TRUE)
  trna <- subset(all_data, cat == "tRNA")
  trna %>%
    arrange(desc(count)) %>%
    head()
}
