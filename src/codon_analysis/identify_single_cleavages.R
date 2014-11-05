# Identify single site hits
# Some code lifted from diffuse_discrete_scatter.R

# Set working directory and variables
library(dplyr)
setwd("~/projects/5OH/results/20140415_rep1")
table <- "ATTGGC_S4.intersect.tab"
min.count <- 25
num.hits <- 40
output.file <- "top_discrete_cleavages40.tab"

# Read in data and subset on mRNA with ample counts
df <- read.table(table, header=TRUE)
df <- tbl_df(df)
df_toPlot <- subset(df,cat=="mRNA" & count >= min.count & ! gene == "RPL41B")

# Ratio of total counts over number of sites
total_ratio <- df_toPlot %.%
  group_by(gene) %.%
  summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())

# Select genes with single hits
one_hit <- subset(total_ratio, total_counts == ratio)

# Intersect top single cleavage genes with original dataset to get coordinates of cleavage
top_one_hits <- subset(total_ratio, total_counts == ratio) %.%
  arrange(desc(total_counts)) %.%
  head(num.hits)
top_one_hits_sites <- subset(df_toPlot, gene %in% top_one_hits$gene)

# Output top hits as bedfile
df_to_write <- top_one_hits_sites %.% 
  select(chr, start, stop, gene, count) %.% 
  arrange(desc(count))

write.table(df_to_write, file = output.file, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)