library(dplyr)
library(reshape2)
library(ggplot2)

# Plots RNA types in bar form
barplot_RNAtypes <- function(data,sample_name) {
  bar_data <- data %.%
    group_by(gene,cat) %.%
    summarize(num=n())
  
  plot.title="Types of RNA Identified"
  
  gp <- ggplot(data=bar_data, aes(x=cat, fill=factor(cat)))
  gp + geom_bar(position="stack") +
    ggtitle(bquote(atop(bold(.(plot.title)),
                        atop(italic(.(sample_name)), "")))) +
    labs(x="RNA Categories",y="Number of unique RNAs") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
    #scale_y_continuous(breaks=c(2,4,6,8,10))
  #ggsave("WT-DMSO-rep1_RNAtypes.png")
}

piechart_RNAtypes <- function(data,sample_name) {
  
  # Summary by number of sites per gene
  #pie_data <- data %.%
    #group_by(cat,gene) %.%
    #summarize(num=n()) #%.% # first summary allows to count number of sites
    #group_by(cat) %.% # second summarize allows to count number of genes
    #summarize(num=n())
  
  # OR Summary by number of genes per category
  #pie_data <- data %.%
  #  group_by(cat,gene) %.%
  #  summarize(num=n()) %.% 
  #  group_by(cat) %.%
  #  summarize(num=n())
  
  # OR Summary by number of counts per gene
  #pie_data <- data %.%
  #  group_by(gene,cat) %.%
  #  summarize(gene_count=sum(count))
  
  # OR Summary by number of counts per category
  pie_data <- data %.%
    group_by(gene,cat) %.%
    summarize(gene_count=sum(count))
  pie_data <- data.frame(pie_data) %.%
    group_by(cat) %.%
    summarize(cat_counts=sum(gene_count))
  
  plot.title="Types of RNA Identified"
  
  #slices <- pie_data$gene_count
  #slices <- pie_data$num
  slices <- pie_data$cat_counts
  lbls <- pie_data$cat
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct) # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main=plot.title)
  
  #gp <- ggplot(data=pie_data, aes(x=cat,y=num,fill=factor(cat)))
  #gp + geom_bar(position="stack") +
  #  ggtitle(bquote(atop(bold(.(plot.title)),
  #                     atop(italic(.(sample_name)), "")))) +
  #  labs(x="RNA Categories",y="Number of unique RNAs") +
  #  theme(legend.position = "none",
  #        axis.text.x = element_text(angle = 45, hjust = 1))
  #scale_y_continuous(breaks=c(2,4,6,8,10))
  #ggsave("WT-DMSO-rep1_RNAtypes.png")
  
  #pie_data %.% arrange(desc(cat_counts))
}

ggpie_RNAtypes <- function(data,sample_name) {
  
  # Summary by number of sites per gene
  #pie_data <- data %.%
  #group_by(cat,gene) %.%
  #summarize(num=n()) #%.% # first summary allows to count number of sites
  #group_by(cat) %.% # second summarize allows to count number of genes
  #summarize(num=n())
  
  # OR Summary by number of genes per category
  #pie_data <- data %.%
  #  group_by(cat,gene) %.%
  #  summarize(num=n()) %.% 
  #  group_by(cat) %.%
  #  summarize(num=n())
  
  # OR Summary by number of counts per gene
  #pie_data <- data %.%
  #  group_by(gene,cat) %.%
  #  summarize(gene_count=sum(count))
  
  # OR Summary by number of counts per category
  pie_data <- data %.%
    group_by(gene,cat) %.%
    summarize(gene_count=sum(count))
  pie_data <- data.frame(pie_data) %.%
    group_by(cat) %.%
    summarize(cat_counts=sum(gene_count))
  pie_data <- data.frame(pie_data) %.%
    mutate(sample=0)
  
  plot.title="Percentage Of Read Counts by RNA Type"
  
  gp <- ggplot(data=pie_data, aes(x=sample,y=cat_counts,fill=factor(cat)))
  gp + geom_bar(width=1,stat="identity",position="stack") +
    coord_polar(theta="y") +
    xlab('') +
    ylab('') +
    labs(fill='RNA Type') +
  ggtitle(bquote(atop(bold(.(plot.title)),
                       atop(italic(.(sample_name)), "")))) #+
  #  labs(x="RNA Categories",y="Number of unique RNAs") +
  #  theme(legend.position = "none",
  #        axis.text.x = element_text(angle = 45, hjust = 1))
  #scale_y_continuous(breaks=c(2,4,6,8,10))
  #ggsave("WT-DMSO-rep1_RNAtypes.png")
  
  #pie_data %.% arrange(desc(cat_counts))
}


wt <- read.delim("ATTGGC_S4.intersect.tab", header=TRUE)
wt_sap <- read.delim("GTAGCC_S8.intersect.tab", header=TRUE)

highly_abundant <- c("gene","CDS_Verified","noncoding_exon",#"mRNA","rRNA",
                     "CDS_Dubious","CDS_Uncharacterized","CDS_Verified%7Csilenced_gene")
less_interesting <- c("Y_prime_element","internal_transcribed_spacer_region","external_transcribed_spacer_region")
other_annotations <- c("non_transcribed_region","pseudogene",
                       "five_prime_UTR_intron","ncRNA","intron","snoRNA")

bar_genes <- subset(wt,
                    ! cat %in% highly_abundant & count > 9 & ! cat %in% less_interesting)

#barplot_RNAtypes(bar_genes,"WT DMSO rep1")
piechart_RNAtypes(filter(bar_genes,! cat %in% other_annotations),"WT DMSO rep1")

ggpie_RNAtypes(subset(wt,
                      count > 9 &
                        ! cat %in% highly_abundant &
                        ! cat %in% less_interesting &
                        ! cat %in% other_annotations)
               ,"WT DMSO rep1")
###