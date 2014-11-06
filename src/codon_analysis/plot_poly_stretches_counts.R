# Barplot of poly-aa count data
# refusing to dodge.  screw you R.
setwd("~/devel/5OH/src/codon_analysis")

library(dplyr)
library(Cairo)

data <- read.table("poly_stretches_counts.tab", header=TRUE)
data <- subset(data, num.rep > 5)

gp <- ggplot(data, aes(mismatch, num.genes, fill = num.rep)) + 
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~polyaa) +
  theme_bw()

gp_single <- ggplot(subset(data, ! polyaa %in% c("[DE]", "[RK]")),
                    aes(mismatch, num.genes, fill = num.rep)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~polyaa) + 
  theme_bw()

ggsave(filename = "polyaa_counts.pdf", 
       plot = gp,
       device = CairoPDF)

ggsave(filename = "polyaa_counts_singleaa.pdf", 
       plot = gp_single,
       device = CairoPDF)