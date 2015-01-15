# graph_PARS_vs_CPMs.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH-seq Methods Paper
# 5OH pipeline - Plotting PARS scores vs. CPMs

library(ggplot2)
library(dplyr)

# legacy code for low-throughput version...
#datadir <- "~/projects/5OH/results/methodspaper/pars"
#setwd(datadir)
#sampletable <- "SP8.assembly-sacCer1.align-uniq.strand.all.CPMs-PARS.tab"

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("usage: Rscript graph_PARS_vs_CPMs.R sampletable")
}
sampletable = args[1]

sampletable = "/vol2/home/speach/projects/5OH/results/methodspaper/pars/SP8.assembly-sacCer1.align-uniq.strand.all.CPMs-PARS.tab"

# remove .tab extension
samplename <- substr(sampletable, 1, nchar(sampletable) - 4) 
data <- read.table(sampletable, header=TRUE)
print(head(data,10))

for (cutoff in c(1,10,25,100)) {
  toPlot <- subset(data, CPMs >= cutoff)
  corrs <- cor.test(toPlot$CPMs, toPlot$PARS, method = "pearson")
  
  # calculate correlations
  ggplot(toPlot,aes(x=CPMs, y=PARS)) +
    geom_point() +
    scale_x_log10() +
    theme_bw() +
    geom_text(data=NULL,
              aes(label = paste("min.CPMs=", cutoff, "\n",
                                "n.sites=", corrs$parameter, "\n", 
                                "pearson.R=", signif(corrs$estimate,2), "\n", 
                                "p.value=", signif(corrs$p.value,2), 
                                sep="")),
              x=Inf, y=Inf, hjust=1, vjust=1)
  
  outputpdf <- paste(samplename, ".cutoff-", cutoff, ".pdf", sep="")
  print(outputpdf)
  ggsave(filename=outputpdf,
         useDingbats=FALSE)
}
