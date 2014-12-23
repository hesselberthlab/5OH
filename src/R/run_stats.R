# Analyzes basic sequencing run stats
# - Comparison of Aligned Reads and UMI correction
# - Outputs averages of alignment, umireads, and sites

library(ggplot2)
library(dplyr)
library(reshape2)

genomesize=12157105  #yeast; length of genome mapped to.

stats <- read.table("alignment_stats.txt",header=TRUE)
stats$sample <- c("S1","S2","S3","S4","S5","S6","S7","S8")
stats <- mutate(stats,aligned.pct=aligned/processed,umi.pct=umireads/aligned,site.pct=sites/genomesize)
attach(stats)

# Compare experimental subpopulations
sap_samples <- c("S5","S6","S7","S8")
sap <- subset(stats,sample %in% sap_samples)
no_sap <- subset(stats, ! sample %in% sap_samples)

t.test(sap$aligned,no_sap$aligned) # significantly different
t.test(sap$umireads,no_sap$umireads) # not significantly different

# print average alignment, umireads, sites for given experimental set
output_stats <- function(stats) {
  print ("Mean Aligned Reads:")
  print (mean(stats$aligned))
  print ("Mean UMI Reads:")
  print (mean(stats$umireads))
  print ("Unique Mapping sites")
  print (mean(stats$sites))
  print ("Percentage Aligned Reads:")
  print (mean(stats$aligned.pct))
  print ("Percentage of UMI reads from total")
  print (mean(stats$umi.pct))
  print ("Percentage coverage of genome")
  print (mean(stats$site.pct))
}

output_stats(stats)

# Generate graph of mapping for samples
to_plot <- data.frame(x=sample,Aligned.Reads=aligned,UMI.Reads=umireads)#,Sites=sites)
melted<-melt(to_plot, id="x")

ggplot(melted,aes(x=x,y=value,fill=variable)) + 
  geom_bar(stat="identity",position = "dodge", alpha=.3) +
  ggtitle("UMI correction of Aligned Reads") +
  xlab("Samples") +
  ylab("Reads") +
  theme(legend.title = element_blank())
ggsave("UmiStats.png")