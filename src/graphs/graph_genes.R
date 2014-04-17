# Create UCSC-like graph
# Inputs: Gene coordinates

library(ggplot2)
library(dplyr)

wt <- read.table("ATTGGC_S4.intersect.tab",header=TRUE)
wt_mRNA <- subset(wt,cat=="mRNA")
tm <- read.table("CACTGT_S3.intersect.tab",header=TRUE)
tm_mRNA <- subset(tm,cat=="mRNA")

graph_gene <- function(df,sample_name,gene.data) {
  chrom <- gene.data[1]
  gene.start <- as.numeric(gene.data[2])
  gene.stop <- as.numeric(gene.data[3])
  gene.name <- gene.data[4]  
  df <- subset(df,chr==chrom & start >= gene.start & stop <= gene.stop)
  maxcount <- max(df$count)

  #df2 <- subset(df2,chr==chrom & start >= gene.start & stop <= gene.stop)
  
  ggplot(df,aes(x=start,y=count)) +
    geom_bar(stat="identity") +
    #geom_bar(data=df2,stat="identity") +
    #ggtitle(bquote(atop(bold(.(gene.name)),
    #                    atop(italic(.(sample_name)), "")))) + 
    #labs(x=paste("Genomic Coords" , " (", chrom, ")" , sep=""),
    #     y="Counts Per Million (CPM)") +
    labs(x="",y="") +
    annotate("text", x = gene.stop-50, y=maxcount-10, label = sample_name, size=3) +
    annotate("text", x = gene.start, y=maxcount-10, label = gene.name, size=3) +
    annotate("text", x = gene.start, y=0, label = chrom, size=3) +
    xlim(gene.start, gene.stop) +
    theme_bw()
  
  filename = paste((paste(sample_name,gene.name,sep="_")),"png",sep=".")
  #BIG
  ggsave(filename,height=2,width=6)
  #small
  #ggsave(filename,height=2,width=3)
}

ADE8 = c("chr4",1288215,1288859,"ADE8","-")
graph_gene(wt_mRNA,"WT-DMSO-rep1",ADE8)
graph_gene(tm_mRNA,"WT-Tm-rep1",ADE8)

RPS31 = c("chr12",498947,499405,"RPS31","+")
graph_gene(wt_mRNA,"WT-DMSO-rep1",RPS31)
graph_gene(tm_mRNA,"WT-Tm-rep1",RPS31)

MDH1 = c("chr11",279123,280127,"MDH1","+")
graph_gene(wt_mRNA,"WT-DMSO-rep1",MDH1)

TDH1 = c("chr10", 338271, 339269, "TDH1",	"+")
graph_gene(wt_mRNA,"WT-DMSO-rep1",TDH1)