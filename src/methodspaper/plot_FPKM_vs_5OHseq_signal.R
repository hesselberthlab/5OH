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

datadir <- "~/projects/5OH/results/methodspaper"
setwd(datadir)

samplename <- "SP8.assembly-sacCer1.align-all"
sampletable <- paste(samplename, ".FPKM_sum.tab", sep="")


df <- read.table(sampletable, header=FALSE,
                 col.names = c("chrom", "start", "stop", "gene", "FPKM", "strand", "OH.signal"))

highFPKMgenes <- subset(df, FPKM > 16000)
highFPKMgenes %>% arrange(desc(FPKM))

high5OHgenes <- subset(df, OH.signal > 3000 & FPKM < 20000)
high5OHgenes %>% arrange(desc(OH.signal))

p <- ggplot(subset(df, OH.signal != 0 & FPKM != 0), aes(x=OH.signal, y=FPKM)) +
  geom_point() +
  stat_smooth(method=lm, se=FALSE) +
  geom_text(data=highFPKMgenes, mapping=aes(label=gene),
            size=3.5,hjust=-0.1,vjust=0) +
  geom_text(data=high5OHgenes, mapping=aes(label=gene),
            size=3.5,hjust=0.5,vjust=-0.5) +
  xlab("5OH-Seq signal (CPMs)") +
  #scale_x_log10() +
  #scale_y_log10() + 
  theme_bw()

df.lm = lm(FPKM ~ OH.signal, data=df)
summary(df.lm)