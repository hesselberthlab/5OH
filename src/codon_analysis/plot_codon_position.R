library(ggplot2)
library(Cairo)

setwd("~/devel//5OH/src//codon_analysis")
pdf.filename <- "codon_position_top40singlesitegenes.pdf"
pdf.facet.filename <- "codon_position_top40singlesitegenes_strandfacet.pdf"

data <- read.table("codons40.bed", col.names = c("chrom","start","stop","gene","position","strand"))
gp <- ggplot(data, aes(factor(position))) + 
  geom_bar() +
  theme_bw() +
  xlab("Codon Position") +
  ylab("Number of genes")

ggsave(filename = pdf.filename, 
       plot = gp,
       device = CairoPDF)

ggsave(filename = pdf.facet.filename, 
       plot = gp + facet_wrap(~strand),
       device = CairoPDF)

