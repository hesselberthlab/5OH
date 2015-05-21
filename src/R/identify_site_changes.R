# identify_site_changes.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH pipeline - Plotting site changes between two datasets
# Plots counts of individual 5OH sites in two datasets on opposing axes
# Colors outliers above a fold-change threshold

library(dplyr)
library(ggplot2)
library(Cairo)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
  stop("usage: Rscript identify_site_changes.R sample1.intersect.tab sample1.name sample2.intersect.tab sample2.name min.count fold.change output.dir")
}

# Parse command line arguments; end process if there is an error
x.table = args[1]
x.name = args[2]
y.table = args[3]
y.name = args[4]
min.count = args[5]
fold.change = as.integer(args[6])
output.dir = args[7]

dfx <- read.table(x.table,header=TRUE) %.%
  subset(count >= min.count)
dfy <- read.table(y.table,header=TRUE)  %.%
  subset(count >= min.count)


# Fitler dataframe by mRNA category
get_mRNA <- function(df) {
  df <- df %.%
    filter(cat=="mRNA") %.%
    select(gene,site,count)
  
  return(df)
}

# Generate plot of sites in df.y vs. df.x
plot_site_counts <- function(df,x_changes,y_changes,dfx.label,dfy.label) {
  gp <- ggplot(df, aes(x=count.x,y=count.y)) +
    geom_point() +
    xlim(0,500) + ylim(0,500) +
    geom_smooth(method = "lm", color="black", formula = y ~ x) +
    geom_point(data=x_changes,color="red") +
    geom_point(data=y_changes,color="darkgreen") +
    xlab(dfx.label) + ylab(dfy.label) + 
    geom_text(data=filter(x_changes,count.x > 150),mapping=aes(label=gene),size=4) + 
    geom_text(data=filter(y_changes,count.y >= count.x*10,count.y>100),
              mapping=aes(label=gene),size=4,
              vjust=1)
  
  return (gp)
}

# Grab statistics on comparison
# TODO: Change to output somewhere in a useful format
stats_for_site_comparison <- function(x,y,dfx.label,dfy.label) {
  x_sites <- x %.% 
    filter(count.x > 39) %.% 
    group_by(gene) %.% 
    arrange(desc(count.x/count.y))
  x_genes <- x_sites %.%
    summarize(n())
  
  y_sites <- y %.% 
    filter(count.y > 39) %.% 
    group_by(gene) %.% 
    arrange(desc(count.y/count.x))
  y_genes <- y_sites %.%
    summarize(n())
  
  print (paste("Num sites increased in ",dfx.label,": ",nrow(x_sites),sep=""))
  print (paste("Num genes increased in ",dfx.label,": ",nrow(x_genes),sep=""))
  
  print (paste("Num genes increased in ",dfy.label,": ",nrow(y_genes),sep=""))
  print (paste("Num sites increased in ",dfy.label,": ",nrow(y_sites),sep=""))
  
  return (list(x=x_sites,y=y_sites))
}

# Merge dataframes and pass to plotting function
site_plot <- function(dfx,dfy,dfx.label,dfy.label,fold) {
  dfx_mRNA <- get_mRNA(dfx)
  dfy_mRNA <- get_mRNA(dfy)
  
  dfs_join <- inner_join(dfx_mRNA,dfy_mRNA, by=c("gene","site"))
  x_change <- subset(dfs_join, count.x >= count.y * fold)
  y_change <- subset(dfs_join, count.y >= count.x * fold)
  
  plot_site_counts(dfs_join,x_change,y_change,dfx.label,dfy.label)
}

# Merge dataframes and pass to stats function
site_stats <- function(dfx,dfy,dfx.label,dfy.label,fold) {
  dfx_mRNA <- get_mRNA(dfx)
  dfy_mRNA <- get_mRNA(dfy)
  
  dfs_join <- inner_join(dfx_mRNA,dfy_mRNA, by=c("gene","site"))
  x_change <- subset(dfs_join, count.x >= count.y * fold)
  y_change <- subset(dfs_join, count.y >= count.x * fold)
  
  stats_for_site_comparison(x_change,y_change,dfx.label,dfy.label)
}

# Plot data and save file
isc_gp <- site_plot(dfx,dfy,x.name,y.name,fold.change)
site_stats(dfx,dfy,x.name,y.name,fold.change)

sample.name <- paste(x.name, ".vs.", y.name, sep='')

pdf.filename <- paste(output.dir, '/', '5OH.sitechanges.',
                      min.count, "-mincount.",
                      fold.change, "-foldchange.",
                      sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, 
       plot = isc_gp,
       device = CairoPDF)
