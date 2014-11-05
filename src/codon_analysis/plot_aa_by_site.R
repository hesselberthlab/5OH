library(ggplot2)
library(Cairo)

setwd("~/devel//5OH/src//codon_analysis")
pdf.filename <- "aa_site_top20ss.pdf"

data <- read.table("top20ss_aa.tab", fill = TRUE,
                   col.names = c("position","aa","count"))
data <- data %.% arrange(aa) %.% arrange(position)
black <- "#000000"
red <- "#FF0000"
blue <- "#0000FF"
gp <- ggplot(data, aes(x=position,y=count, fill=aa)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c(black, black, black, red, red,
                             black, black, blue, black, blue,
                             black, black, black, black, black,
                             blue, black, black, black, black, black, black, black)) +
  xlab("Position") +
  ylab("Number of amino acids")

ggsave(filename = pdf.filename, 
       plot = gp,
       device = CairoPDF)

