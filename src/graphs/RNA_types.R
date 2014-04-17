# Plots RNA types in bar form
barplot_RNAtypes <- function(data,sample_name) {
  bar_data <- data %.%
    group_by(gene,cat) %.%
    summarize(num=n())
  
  plot.title="Types of RNA Identified"
  
  gp <- ggplot(data=bar_data, aes(x=cat,fill=factor(cat)))
  gp + geom_bar(position="stack") +
    ggtitle(bquote(atop(bold(.(plot.title)),
                        atop(italic(.(sample_name)), "")))) +
    labs(x="RNA Categories",y="Number of unique RNAs") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks=c(2,4,6,8,10))
  ggsave("WT-DMSO-rep1_RNAtypes.png")
}

wt <- read.delim("ATTGGC_S4.intersect.tab", header=TRUE)
wt <- subset(wt, count > 4)

highly_abundant <- c("gene","CDS_Verified","noncoding_exon","mRNA","rRNA",
                     "CDS_Dubious","CDS_Uncharacterized","CDS_Verified%7Csilenced_gene")
less_interesting <- c("Y_prime_element","internal_transcribed_spacer_region","external_transcribed_spacer_region")

bar_genes <- subset(wt,
                    ! cat %in% highly_abundant & count > 9 & ! cat %in% less_interesting)
barplot_RNAtypes(bar_genes,"WT DMSO rep1")
###