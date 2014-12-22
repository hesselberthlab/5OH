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

# Directory containing results
setwd("~/projects/5OH/results/methodspaper/")

# Useful library prefixes
sample_libraries = c("SP5", "SP6", "SP7", "SP8", "SP9", "SP10", "SP11", "SP12", "SP13", "SP14", "SP15", "SP16", "SP17", "SP18", "SP19", "SP27")
noSAP_libraries = c("SP5", "SP6", "SP7", "SP8", "SP13", "SP14", "SP15", "SP27")
SAP_libraries = c("SP9", "SP10", "SP11", "SP12", "SP16", "SP17", "SP18", "SP19")
control_libraries = c("SP20", "SP28", "SP29")
all_libraries = c("SP5", "SP6", "SP7", "SP8", "SP9", "SP10", "SP11", "SP12", "SP13", "SP14", "SP15", "SP16", "SP17", "SP18", "SP19", "SP27", "SP20", "SP28", "SP29")

# Alignment specifics
assembly <- "assembly-sacCer1"
alignment <- "align-all"

# Sample specifics
workingdir = "~/projects/5OH/results/methodspaper/"
sample.name = "SP8"
sample.descr = "WT_DMSO_Rep1"
output.dir = paste(workingdir, "rgraphs", sep="")
fileprefix <- paste(sample.name, assembly, alignment, sep=".")
table <- paste("/rgraphs/", fileprefix, ".gff.intersect", sep="")


# Function to calculate single cleavage sites
get_one_hit_genes <- function(df){
  total_ratio <- df %>%
    group_by(gene) %>%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n()) %>%
    arrange(desc(total_counts))
  
  one_hit_genes <- subset(total_ratio, total_counts == ratio)
  
  return (one_hit_genes)
}

# Function to generate diffuse_discrete_scatterplot
diffuse_plot <- function(df,sample.descr,lines=TRUE){
  plot.title = "Diffuse and Discrete Cleavage"
  
  total_ratio <- df %>%
    group_by(gene) %>%
    summarize(total_counts=sum(count), ratio=total_counts/n(), sites=n())
  
  maxy = max(total_ratio$ratio)
  minx = min(total_ratio$total_counts)
  
  one_hit <- subset(total_ratio, total_counts == ratio)
  top_one_hits <- subset(total_ratio, total_counts == ratio) %>%
    arrange(desc(total_counts)) %>%
    head(10)
  print(one_hit)
  
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

setwd(workingdir)

# Parse input file, create subset above minimum count
min.count = 25
df <- read.table(table, header=TRUE)
df <- tbl_df(df)
df_toPlot <- subset(df,cat=="mRNA" & count >= min.count & ! gene == "RPL41B")

# Plot data and save file
dd_gp <- diffuse_plot(df_toPlot, sample.descr)
dd_gp
one_hits <- get_one_hit_genes(df_toPlot)
one_hits

pdf.filename <- paste(output.dir, '/', '5OH.ddscatter.',
                      min.count, "-mincount.",
                      sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, 
       plot = dd_gp,
       device = CairoPDF)
