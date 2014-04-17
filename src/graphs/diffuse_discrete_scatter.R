# Look at Diffuse and discrete Cleavages

library(ggplot2)
library(dplyr)
library(reshape2)

# A plot to look at diffuse and discrete cleavages
# Discrete appears in top left (ie, RPS31)
# Diffuse appears at bottom right (ie, CBF5)
plot_diffuse <- function(data,sample_name){
  plot.title = "Diffuse and Discrete Cleavage"
  
  total_ratio <- data %.%
    group_by(gene) %.%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())
  
  maxy = max(total_ratio$ratio)
  minx = min(total_ratio$total_counts)
  
  #tr_labels <- subset(total_ratio, total_counts>200 | ratio > 100)
  #tr_labels <- subset(total_ratio, ratio > 70)
  #tr_labels <- total_ratio
  tr_labels <- subset(total_ratio, gene %in% c(#"RPS31","ADE8","TPI1","YNK1","TPI1",
                                               "MDH1","TDH1","RPL40A"))
  one_hit <- subset(total_ratio, total_counts == ratio)
  top_one_hits <- subset(total_ratio, total_counts == ratio) %.%
    arrange(desc(total_counts)) %.%
    head(5)
  #tr_labels <- subset(one_hit, total_counts > 64)
  
  gp <- ggplot((total_ratio), aes(y=ratio, x=total_counts))
  gp + geom_point(size=2.5) +
    geom_text(data=tr_labels, mapping=aes(label=gene),
              size=3.5,hjust=1.1, vjust=0) +
    #geom_text(data=top_one_hits, mapping=aes(label=gene),
    #          size=3.5,hjust=1,vjust=0) +
    ggtitle(bquote(atop(bold(.(plot.title))))) +
    labs(x="Total Number of UMI-corrected reads",y="Ratio: Total Reads / Number of sites") +
    scale_x_log10() +
    annotate("text", x = minx+10, y=maxy, label = sample_name, size=4) +
    #geom_line(data=one_hit, color="red") +
    #geom_line(data=subset(total_ratio,total_counts/2==ratio), color="red") +
    #geom_line(data=subset(total_ratio,total_counts/3==ratio),color="red") +
    #geom_line(data=subset(total_ratio,total_counts/4==ratio),color="red") +
    #geom_line(data=subset(total_ratio,total_counts/5==ratio),color="red") +
    #geom_line(data=subset(total_ratio,total_counts/6==ratio),color="red") +
    theme_bw()
}
tm <- read.table("CACTGT_S3.intersect.tab",header=TRUE)
wt <- read.table("ATTGGC_S4.intersect.tab",header=TRUE)
plot_diffuse(subset(tm,cat=="mRNA" & count > 19 & ! gene == "RPL41B"), "WT +Tm; Rep1")
plot_diffuse(subset(wt,cat=="mRNA" & count > 24 & ! gene == "RPL41B"),"WT-DMSO-Rep1")



### 
plot_diffuse_lee_perfect <- function(data,sample_name){
  total_ratio <- data %.%
    group_by(gene) %.%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())
  
  #tr_labels <- subset(total_ratio, total_counts>1000 | ratio > 75)
  tr_labels <- subset(total_ratio, ratio > 70)
  #tr_labels <- total_ratio
  one_hit <- subset(total_ratio, total_counts == ratio)
  top_one_hits <- subset(total_ratio, total_counts == ratio) %.%
    arrange(desc(total_counts)) %.%
    head(5)
  #tr_labels <- subset(one_hit, total_counts > 64)
  
  gp <- ggplot((total_ratio), aes(y=ratio, x=total_counts))
  gp + geom_point() +
    #geom_text(data=tr_labels, mapping=aes(label=gene),
    #          size=3.5,hjust=-.1, vjust=0) +
    geom_text(data=top_one_hits, mapping=aes(label=gene),
              size=3.5,hjust=1,vjust=0) +
    ggtitle(expression(atop(bold("Diffuse and Discrete Cleavage"),
                            atop(italic("WT DMSO: 20140107 S4"), "")))) +
    labs(x="Total Number of UMI-corrected reads, log(10)",y="Ratio: Total Reads / Number of sites") +
    scale_x_log10() +
    geom_line(data=one_hit, color="red") +
    geom_line(data=subset(total_ratio,total_counts/2==ratio), color="red") +
    geom_line(data=subset(total_ratio,total_counts/3==ratio),color="red") +
    geom_line(data=subset(total_ratio,total_counts/4==ratio),color="red") +
    geom_line(data=subset(total_ratio,total_counts/5==ratio),color="red") +
    geom_line(data=subset(total_ratio,total_counts/6==ratio),color="red") +
    theme_bw()
}
plot_diffuse_lee_perfect(subset(lee,count>29 & gene!="RPL41B"))
####