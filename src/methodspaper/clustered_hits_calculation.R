# Takes tab-delimitted output file of macs_intersect.sh and counts number of unique mRNA species over some threshold
library(dplyr)

workingdir = "~/projects/5OH/results/methodspaper/"
setwd(workingdir)

# Sample data
sample <- "SP8"
assembly <- "assembly-sacCer1"
alignment <- "align-uniq"
fileprefix <- paste(sample.name, assembly, alignment, sep=".")
filename <- paste("peaks/", fileprefix, "_peaks.intersect.tab", sep="")

df <- read.table(filename, col.names = c("peak.chrom", "peak.start", "peak.stop", "gene", "peak.reads", "site.proportion", "site.chrom", "site.position", "strand"))

df <- df %>%
  mutate(peak.length = peak.stop - peak.start) %>%
  mutate(site = paste(site.chrom, site.position, sep=":")) %>%
  select(gene, site, site.proportion, peak.length, peak.reads) %>%
  arrange(desc(peak.length))

# Removing duplciate entries
df <- df[!duplicated(df$site),]

# List of discrete hits from discrete_hits_indentification.R; want to remove from list of clustered hits to ensure
# no overlap. Currently hardcoded and put together outside of script due to initial manual validation in UCSC browser.
peaklist <- levels(read.table("peaks/peaklist.txt", col.names = c("gene"))$gene)

# Site proportion threshold site to <= 0.15
# Want the highest signal in the gene body to comprise a low proportion of signal due to plethora of other cleavages
subset(df, ! gene %in% peaklist & site.proportion <= 0.15) %>%
  summarize(mean(peak.length)) # Reports average peak length for clustered hits

# Report number of mRNA species with clustered cleavages for this criteria
nrow(subset(df, ! gene %in% peaklist & site.proportion <= 0.15)) 
