# Create UCSC-like graph using merged data frames of samples and
# coordinates of gene to graph.  Supports introns.
# Currently strandedness doesn't work.  :/

library(dplyr)
library(ggplot2)

graph_gene <- function(df,gene.data,intron=FALSE) {
  chrom <- gene.data[1]
  gene.start <- as.numeric(gene.data[2])
  gene.stop <- as.numeric(gene.data[3])
  gene.name <- gene.data[4]
  gene.strand <- gene.data[5]
  
  # Change category call for rRNA vs. mRNA vs. other graphings
  df <- subset(df,chr==chrom & cat=="rRNA" & start >= gene.start & stop <= gene.stop)
  maxcount <- max(df$count)
  base=maxcount-(maxcount/10)

  gp <- ggplot(df,aes(x=start,y=count))
  
  if (gene.strand == "-") {
    start = gene.start
    gene.start = gene.stop
    gene.stop = start
    gp <- gp + scale_x_reverse(limits=c(gene.start-25,gene.stop+25))
  } 
  else {
    gp <- gp + xlim(gene.start-25, gene.stop+25)
  }
  
  fiveprime = gene.start
  threeprime = gene.stop
  
  #df2 <- subset(df2,chr==chrom & start >= gene.start & stop <= gene.stop)
  gp <- gp +
    geom_bar(stat="identity") + facet_grid(sample ~ .) +
    labs(x=paste("Genomic Coords" , " (", chrom, ")" , sep=""),
         y="Counts Per Million (CPM)",
          title=gene.name) +
    theme(axis.title.x=element_text(size=12),
          strip.text.x = element_text(size = 12, colour = "red", angle = 0)) +
    theme_bw() + 
    geom_rect(xmin=fiveprime,
              xmax=threeprime,
              ymin=base,
              ymax=Inf,
              fill="blue",group=NULL,alpha=0.5,data=subset(df,sample=="WT_DMSO_rep1")) +
    geom_text(data=subset(df,sample=="WT_DMSO_rep1"),
                          x = fiveprime - 10,
                          y=maxcount-(maxcount/15),
                          label = "5'", size=5) + 
    geom_text(data=subset(df,sample=="WT_DMSO_rep1"),
              x = threeprime + 10,
              y=maxcount-(maxcount/15),
              label = "3'", size=5)
  
  if (length(intron)==5) {
    intron.start = as.numeric(intron[2])
    intron.stop = as.numeric(intron[3])

    gp <- gp + 
      geom_rect(xmin=intron.start,
                xmax=intron.stop,
                ymin=base,
                ymax=Inf,
                fill="grey",group=NULL,data=subset(df,sample=="WT_DMSO_rep1"))
  }
  
  gp
  
  
  #filename = paste((paste(sample_name,gene.name,sep="_")),"png",sep=".")
  #BIG
  #ggsave(filename,height=2,width=6)
  #small
  #ggsave(filename,height=2,width=3)
}

# have been using all_samples from boxplot code merges

graph_gene(all_samples,HAC1,intron=HAC1_intron)
graph_gene(all_samples,ADE8)
graph_gene(subset(all_samples,strain=="WT"),KAR2)

# Genes I've graphed
ADE8 = c("chr4",1288215,1288859,"ADE8","-")
KAR2 = c("chr10",381327,383375,"KAR2","+")
RDN25 = c("chr12",451786,455181,"RDN25-1","-")
RPL40A = c("chr9",68708,69528,"RPL40A","+")
HAC1 = c("chr6",75179,76147,"HAC1",  "+")
HAC1_intron = c("chr6",75840,76091,"HAC1_intron","+")
TDH1 = c("chr10", 338271, 339269, "TDH1", "+")
MDH1 = c("chr11",279123,280127,"MDH1","+")
RPS31 = c("chr12",498947,499405,"RPS31","+")
RDN25_5prime = c("chr12",451786,455181,"RDN25-1","-")
M19229_28s = c("chr12",454814,455180,"RDN25-1","-")
M19229_28s_5prime = c("chr12",455179,455180,"RDN25-1, 5'end","-")


# Have to tweak graph gene for rRNA
graph_gene(all_samples,M19229_28s_5prime)
