# Plotting alignment results by sample, assembly, alignment method

library(ggplot2)
library(reshape2)
library(dplyr)

setwd("~/projects/5OH/results/methodspaper")
filename <- "alignment_stats.tab"
data <- read.table(filename, header=TRUE) %>%
  select(sample, assembly, alignment, processed, aligned, failed) %>%
  melt(id.vars=c("sample", "assembly", "alignment"))

sample_libraries = c("SP5", "SP6", "SP7", "SP8", "SP9", "SP10", "SP11", "SP12", "SP13", "SP14", "SP15", "SP16", "SP17", "SP18", "SP19", "SP27")
SAP_libraries = c("SP9", "SP10", "SP11", "SP12", "SP16", "SP17", "SP18", "SP19")
control_libraries = c("SP20", "SP28", "SP29")

# Subset data on "true" samples, bowtie -all alignment, sacCer1 genome
samples <- subset(data, sample %in% sample_libraries & alignment == "all" & assembly == "sacCer1")

# Output averages; these values used in paper
samples %>%
  group_by(variable) %>%
  summarize(avg = mean(value))

# Can visualize quick summary plots of alignment data; SP5 is a subpar library
ggplot(samples, aes(x=variable, y=value, fill=sample)) +
  geom_bar(stat="identity", position="dodge")
