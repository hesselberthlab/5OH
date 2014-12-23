# I am actively and frequently editing this file to visualize and interpret data.
# There's code all over the place.  Still... version contol.

library(ggplot2)
library(dplyr)
library(reshape2)

lee <- read.delim("S4.lee.tab", header=TRUE)
all_genes2 <- read.delim("S4.intersect_2.tab", header=TRUE)
all_genes2 <- subset(all_genes2, count > 4)
#genes <- subset(all_genes, count > 19)
#highest <- subset(all_genes, count > 500)

genes2 <- subset(all_genes2,count > 19 & cat != "gene" & cat != "CDS_Verified")
highest2 <- subset(genes2, count >500)

# Code tidbits
#gp <- ggplot(genes2, aes(x=cat, y=count))
#gp + geom_boxplot() +
#  geom_text(data=highest2, mapping=aes(label=gene),
#            size=3.5,hjust=-0.1, vjust=0)

# Bar graph of less abundant RNA types
highly_abundant <- c("gene","CDS_Verified","mRNA","noncoding_exon","rRNA",
                     "CDS_Dubious","CDS_Uncharacterized","CDS_Verified%7Csilenced_gene")
less_interesting <- c("Y_prime_element","internal_transcribed_spacer_region","external_transcribed_spacer_region")

barplot_RNAtypes <- function(data) {
  pie_data <- data %.%
    group_by(gene,cat) %.%
    summarize(num=n())
  
  gp <- ggplot(data=pie_data, aes(x=cat,fill=factor(cat)))
  gp + geom_bar(position="stack") +
    ggtitle(expression(atop(bold("Types of RNA Identified"),
                            atop(italic("WT DMSO: 20140107 S4 pos"), "")))) +
    labs(x="RNA Categories",y="Number of unique RNAs") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

### DO NOT CHANGE
barplot_RNAtypes <- function(data) {
  bar_data <- data %.%
    group_by(gene,cat) %.%
    summarize(num=n())
  
  gp <- ggplot(data=bar_data, aes(x=cat,fill=factor(cat)))
  gp + geom_bar(position="stack") +
    ggtitle(expression(atop(bold("Types of RNA Identified"),
                            atop(italic("WT DMSO: 20140107 S4 pos"), "")))) +
    labs(x="RNA Categories",y="Number of unique RNAs") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}
highly_abundant <- c("gene","CDS_Verified","noncoding_exon","mRNA","rRNA",
                     "CDS_Dubious","CDS_Uncharacterized","CDS_Verified%7Csilenced_gene")
less_interesting <- c("Y_prime_element","internal_transcribed_spacer_region","external_transcribed_spacer_region")
bar_genes <- subset(all_genes2,
                    ! cat %in% highly_abundant & count > 9 & ! cat %in% less_interesting)
barplot_RNAtypes(bar_genes)
###

filter(bar_genes,cat=="snRNA")


# A plot to look at diffuse and discrete cleavages
# Discrete appears in top left (ie, RPS31)
# Diffuse appears at bottom right (ie, CBF5)
plot_diffuse <- function(data,sample_name){
  plot.title = "Diffuse and Discrete Cleavage"
  
  total_ratio <- data %.%
    group_by(gene) %.%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())
  
  #tr_labels <- subset(total_ratio, total_counts>1000 | ratio > 75)
  #tr_labels <- subset(total_ratio, ratio > 70)
  #tr_labels <- total_ratio
  tr_labels <- subset(total_ratio, gene %in% c("RPS31","ADE8","TPI1","TDH1"))
  one_hit <- subset(total_ratio, total_counts == ratio)
  top_one_hits <- subset(total_ratio, total_counts == ratio) %.%
    arrange(desc(total_counts)) %.%
    head(5)
  #tr_labels <- subset(one_hit, total_counts > 64)
  
  gp <- ggplot((total_ratio), aes(y=ratio, x=total_counts))
  gp + geom_point() +
    geom_text(data=tr_labels, mapping=aes(label=gene),
              size=3.5,hjust=-.1, vjust=0) +
    #geom_text(data=top_one_hits, mapping=aes(label=gene),
    #          size=3.5,hjust=1,vjust=0) +
    ggtitle(bquote(atop(bold(.(plot.title)),
                            atop(italic(.(sample_name)), "")))) +
    labs(x="Total Number of UMI-corrected reads, log(10)",y="Ratio: Total Reads / Number of sites") +
    scale_x_log10() +
    geom_line(data=one_hit, color="red") +
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
plot_diffuse(subset(wt,cat=="mRNA" & count > 19 & ! gene == "RPL41B"),"WT DMSO; Rep1")

### DO NOT EDIT ###
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

# Get genes with single cleavage site
get_discrete_hits <- function(data,sample_name){
  ratio_data <- data %.%
    group_by(gene) %.%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())
  
  return (subset(ratio_data,sites==1))
}

genes <- subset(all_genes2, count > 9)

#plot_diffuse(genes)
threshold = 19
mRNA <- subset(all_genes2, count > threshold & cat=="mRNA" & gene!="RPL41B")
plot_diffuse(mRNA)
nrow(get_discrete_hits(mRNA))
#plot_diffuse(pie_genes)





# Indentifying top factors for multiple datasets
bg_cutoff = 10
type="rRNA"
genes_wt <- subset(read.delim("S4.intersect_2.tab", header=TRUE),
                   count >= bg_cutoff & cat==type)
genes_wt_tun <- subset(read.delim("S3.intersect.tab", header=TRUE),
                   count >= bg_cutoff & cat==type)
genes_del_xrn1 <- subset(read.delim("S2.intersect.tab",header=TRUE),
                         count >= bg_cutoff & cat==type)

# Find highest hit for given gene; output to tab file
find_and_write_hits <- function(data,file_name) {
  highest_hits <- data %.%
    group_by(gene) %.%
    summarize(highest.count=max(count)) %.%
    arrange(desc(highest.count))
  
  write.table(highest_hits,file=file_name,
              append=FALSE,quote=FALSE,row.names=FALSE,
              sep="\t")
}

find_and_write_hits(genes_wt,"S4.highesthits.tab")


get_top_hits <- function(data,n) {
  top_hits <- data %.%
    group_by(gene) %.%
    summarize(total_counts=sum(count)) %.%
    arrange(desc(total_counts)) %.%
    head(n)
  return (top_hits)
}

wt_top <- get_top_hits(genes_wt,100)
wt_tun_top <- get_top_hits(genes_wt_tun,100)
xrn1_top <- get_top_hits(genes_del_xrn1,100)

wt_diff_genes <- setdiff(wt_top$gene,wt_tun_top$gene)
subset(wt_top,gene %in% wt_diff_genes)

length(intersect(wt_top$gene, xrn1_top$gene))
length(intersect(wt_top$gene, wt_tun_top$gene))

# Ribosomal RNA summing
# Why are 37-1 and 37-2 slightly different?  Am taking the higher of the two


# Comparing normalization schemes
# CPM, counts/totalcounts * 1000000
cpm = c(6.99599,2.044,1.11382,1.29714,0.655112,0.880518,0.697533,1.54515)
to_plot <- data.frame(x=sample,CPM=cpm)#,Sites=sites)
melted<-melt(to_plot, id="x")

ggplot(melted,aes(x=x,y=value,fill=variable)) + 
  geom_bar(stat="identity",position = "dodge", alpha=.3) +
  ggtitle("UMI correction of Aligned Reads") +
  xlab("Samples") +
  ylab("Reads") +
  theme(legend.title = element_blank())
  

# Comparing cleavage sites
tm <- read.table("CACTGT_S3.intersect.tab",header=TRUE)
wt <- read.table("ATTGGC_S4.intersect.tab",header=TRUE)





