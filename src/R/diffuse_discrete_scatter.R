# diffuse_discrete_scatter.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH pipeline - Plotting diffuse and discrete cleavages
# Makes scatterplot of (ratio of total reads / number of sites) vs. (total number of reads) for each gene
# Discrete cleavages appear in top left (ie, RPS31); diffuse cleaveages appears at bottom right (ie, CBF5)

# TODO: Currently, output of graphs is not like that seen in the interactive terminal, even when (presumably!)
# operating on same data.  Figure this shit out.

library(ggplot2)
library(dplyr)
library(Cairo)

# Parse command line arguments; end process if there is an error
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("usage: Rscript diffuse_discrete_scatter.R sample.intersect.tab sample.name sample.descr min.count output.dir")
}

table = args[1]
sample.name = args[2]
sample.descr = args[3]
min.count = args[4]
output.dir = args[5]

# Parse input file, create subset above minimum count
df <- read.table(table, header=TRUE)
df <- tbl_df(df)
df_toPlot <- subset(df,cat=="mRNA" & count >= min.count & ! gene == "RPL41B")

# Function to generate diffuse_discrete_scatterplot
diffuse_plot <- function(df,sample.descr,lines=TRUE){
  plot.title = "Diffuse and Discrete Cleavage"
  
  total_ratio <- df %.%
    group_by(gene) %.%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())
  
  maxy = max(total_ratio$ratio)
  minx = min(total_ratio$total_counts)
  
  one_hit <- subset(total_ratio, total_counts == ratio)
  top_one_hits <- subset(total_ratio, total_counts == ratio) %.%
    arrange(desc(total_counts)) %.%
    head(4)
  
  gp <- ggplot((total_ratio), aes(y=ratio, x=total_counts)) +
    geom_point(size=2.5) +
    geom_text(data=top_one_hits, mapping=aes(label=gene),
              size=3.5,hjust=1,vjust=0) +
    ggtitle(bquote(atop(bold(.(plot.title))))) +
    labs(x="Total Number of UMI-corrected reads",y="Ratio: Total Reads / Number of sites") +
    scale_x_log10() +
    annotate("text", x = minx+10, y=maxy, label = sample.descr, size=4) +
    theme_bw()
  
  if (lines==TRUE) {
    gp <- gp +
      geom_line(data=one_hit, color="red") +
      geom_line(data=subset(total_ratio,total_counts/2==ratio), color="red") +
      geom_line(data=subset(total_ratio,total_counts/3==ratio),color="red") +
      geom_line(data=subset(total_ratio,total_counts/4==ratio),color="red") +
      geom_line(data=subset(total_ratio,total_counts/5==ratio),color="red")
  }
  
  return (gp)
}

# Plot data and save file
dd_gp <- diffuse_plot(df_toPlot, sample.descr)

pdf.filename <- paste(output.dir, '/', '5OH.ddscatter.',
                      min.count, "-mincount.",
                      sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, 
       plot = dd_gp,
       device = CairoPDF)
