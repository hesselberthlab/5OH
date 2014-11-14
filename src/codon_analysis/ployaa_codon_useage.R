library(ggplot2)
library(dplyr)
library(reshape2)
library(Cairo)

####################################
# Codon Useage in polyaa stretches #
####################################

motif_plot <- function(data, motif_name, pdf.filename) {
  data <- subset(data, motif == motif_name)
  long_data <- melt(data, id.vars=c("gene","motif"))
  colnames(long_data)[3] <- "position"
  
  gp <- ggplot(long_data, aes(x=value, fill=position)) +
    geom_bar(position="dodge") +
    xlab("Codon") +
    ylab("Total") +
    ggtitle(motif_name)
  
  ggsave(filename = pdf.filename,
         plot = gp,
         device = CairoPDF)
}

setwd("~/devel//5OH/src//codon_analysis")
data <- read.table("poly-all_gene-codons.tab", header=TRUE)
colnames(data) <- c("gene","motif","1","2","3","4","5")

motif_plot(data, "E5", "polyE_codon_useage.pdf")
motif_plot(data, "D5", "polyD_codon_useage.pdf")
motif_plot(data, "K5", "polyK_codon_useage.pdf")
motif_plot(data, "R5", "polyR_codon_useage.pdf")


############################################################################
# Plotting Codon Useage in polyaa stretches vs. Known Transcriptome Useage #
############################################################################

summary_data <- melt(data, id.vars=c("gene","motif")) %>%
  subset(motif == "D5" & value %in% c("GAT","GAC")) %>%
  select(motif, value) %>%
  group_by(motif, value) %>%
  summarize(count=n()) %>%
  mutate(freq = count / sum(count), data = "Poly-aa D5") %>%
  select(value, freq, data)

known_codon_usage = read.table("yeast_codons.txt",
                               header=FALSE,
                               col.names = c("motif","value","freq")) %>%
  mutate(data = "Global")

gp <- ggplot(rbind(summary_data, known_codon_usage), aes(x=value, y=freq, fill=data)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("Codon") +
  ylab("Frequency") +
  theme(legend.title = element_blank()) +
  theme_bw()

ggsave(filename = "polyE_codon_useage_vs_known.pdf",
       plot = gp,
       device = CairoPDF)



