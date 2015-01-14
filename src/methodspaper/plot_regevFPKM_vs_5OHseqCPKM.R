# plot_regevFPKM_vs_5OHseqFPKM.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH-seq Methods Paper
# Plotting Regev RNA vs. 5OH-seq signals, both as FPKM.

library(ggplot2)
library(dplyr)
library(Cairo)

# Initiate empty variable
cor.df <- data.frame()

# Function to plot/save and calculate pearson correlation
plot_and_save <- function(df, cutoff) {
  # Parse based on cutoff
  df_toPlot <- df %>%
    subset(CPKM >= cutoff) %>%
    group_by(gene) %>%
    summarize(CPKM = max(CPKM), RPKM = max(RPKM))
  
  gp <- ggplot(df_toPlot, aes(x=CPKM, y=RPKM)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    xlab("5OH-Seq (CPKM)") +
    ylab("RNA-Seq (FPKM)")
  
  # Save plot
  ggsave(plot = gp,
         filename = paste("regevFPKM.vs.5OHseqCPKM.",cutoff,"-cutoff.pdf",sep=""),
         height = 5,
         width = 6,
         useDingbats=FALSE)
  
  # Calculate correlations and add to table
  cortest <- cor.test(df_toPlot$CPKM, df_toPlot$RPKM, alternative="two.sided", method="pearson")
  current.cor.df <- data.frame(cutoff, cortest$p.value,cortest$estimate,cortest$parameter)
  
  # Oh this is just the coolest; writing to variable outside of current scope :)
  cor.df <<- rbind(cor.df, current.cor.df)
}

# Read in data
datadir <- "~/projects/5OH/data/methodspaper/sacCer1/expression"
setwd(datadir)
sampletable <- "SP8_CPKMs.vs.regev_FPKMs.tab"
df <- read.table(sampletable, header=TRUE)

# Loop through count cutoffs to generate graphs and save images.
for (cutoff in c(1, 10, 25, 100)) {
  plot_and_save(df, cutoff)
}

# Accessing values of loop, and then handing them off to Prism, cause life.
cor.df