# plot_FPKM_vs_5OHseq_signal.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH-seq Paper
# Plot FPKM from /vol2/home/speach/ref/genomes/sacCer1/exp.fpkm.bed vs. 5OH-seq signal
# Input data is from compare_fpkm_data_and_CPMs.sh; ie, a mapping of 5OH-seq stranded bedgraphs to FPKM BEDfile

library(ggplot2)
library(dplyr)
library(Cairo)

datadir <- "~/projects/5OH/results/methodspaper"
setwd(datadir)

# parse data
samplename <- "SP8.assembly-sacCer1.align-all"
sampletable <- paste(samplename, ".FPKM_max.tab", sep="")
df <- read.table(sampletable, header=FALSE,
                 col.names = c("chrom", "start", "stop", "gene", "FPKM", "strand", "OH.signal"))
filtereddf <- subset(df, FPKM > 0.1 & OH.signal >0) # remove genes without FPKM data

# Interesting things to label
highFPKMgenes <- subset(df, FPKM > 12000 & OH.signal > 10)
highFPKMgenes %>% arrange(desc(FPKM))

high5OHgenes <- subset(df, OH.signal > 1500 & FPKM < 20000)
high5OHgenes %>% arrange(desc(OH.signal))

lowFPKMgenes <- subset(df, FPKM < 0.1)
lowFPKMgenes %>% arrange(desc(FPKM))

discretepeaklist <- levels(read.table("peaks/peaklist.txt", col.names = c("gene"))$gene)
discretehits <- subset(df, gene %in% discretepeaklist)
discrete_figures <- subset(df, gene %in% c("ADE8", "ADK1", "ADH1", "RPS31"))

clusteredpeaklist <- levels(read.table("peaks/clusteredpeaklist.txt", col.names = c("gene"))$gene)
clusteredhits <- subset(df, gene %in% clusteredpeaklist)
clustered_figures <- subset(df, gene %in% c("FUN12", "CBF5"))


full <- ggplot(filtereddf, aes(x=OH.signal, y=FPKM)) +
  geom_point() +
  # Labelled Stuff
  #stat_smooth(method=lm, se=FALSE) +
  #geom_point(data=discretehits, color="blue", size=3) +
  #geom_point(data=clusteredhits, color="red", size=3) +
  #geom_text(data=highFPKMgenes, mapping=aes(label=gene), size=3.5,hjust=-0.1,vjust=0) +
  #geom_text(data=high5OHgenes, mapping=aes(label=gene), size=3.5,hjust=0.5,vjust=-0.5) +
  #geom_text(data=discrete_figures, mapping=aes(label=gene), color="blue",
  #          size=5,hjust=0.5,vjust=-0.5) +
  #geom_text(data=clustered_figures, mapping=aes(label=gene), color="red",
  #          size=5,hjust=0.5,vjust=-0.5) +
  xlab("5OH-Seq signal (CPMs)") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()
full

zoom <- ggplot(subset(filtereddf, FPKM>0.1 & OH.signal>9), aes(x=OH.signal, y=FPKM)) +
  geom_point() +
  xlab("5OH-Seq signal (CPMs)") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()
zoom



df.lm = lm(FPKM ~ OH.signal, data=df)
summary(df.lm)

ggsave(plot = full,
       filename = "SP8.assembly-sacCer1.align-all.FPKM_sum.std.pdf",
       height = 5,
       width = 5,
       useDingbats=FALSE)

ggsave(plot=zoom, filename="SP8.assembly-sacCer1.align-all.FPKM_sum.std.zoom.pdf",
       height=7,
       width=7)