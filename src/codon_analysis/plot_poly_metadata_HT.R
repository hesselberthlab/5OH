# plot_poly_metadata_HT.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH-seq Methods Paper
# Plot poly-aminoacid metadata; high-thruput image generation for sample analysis.

library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(Cairo)

# Parse command line arguments; end process if there is an error
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("usage: Rscript plot_poly_metadata_HT.R input.csv output.pdf")
}
inputfile = args[1]
outputpdf = args[2]

# Full data graphs and subset graphs
read_in_data <- function(filename) {
  data <- read.csv(filename)
  data <- melt(data, id.vars = c("index"),
               variable.name = "polyaa")
  data$value <- as.numeric(data$value)
  
  letters = "([A-Z]|[A-Z][A-Z]|[A-Z][A-Z][A-Z]|[A-Z][A-Z][A-Z][A-Z])"
  numbers = "([0-9]|[0-9][0-9])"
  search = paste("^", letters, numbers, "$", sep="")
  
  parsed_polyaa = as.data.frame(str_match(data$polyaa, search)[,-1])
  colnames(parsed_polyaa) <- c("aa", "num")
  data <- cbind(data, parsed_polyaa)
  
  return(data)
}

plot_polydata <- function(df) {
  ggplot(df, aes(x=index, y=value, color=num)) +
    geom_point() +
    theme_bw() +
    ylab("Average Coverage Per Million Reads") +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75")) +
    xlab("Distance from poly stretch (bp)") +
    facet_wrap(~aa) +
    theme(legend.position="none")
}

plot_polydata_subsets <- function(df,polyaa_list) {
  ggplot(subset(df, polyaa %in% polyaa_list), aes(x=index, y=value, color=polyaa)) +
    geom_point() +
    theme_bw() +
    ylab("Average CPM") +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75")) +
    xlab("Distance from poly stretch (bp)")
}

# It just all looks so cute when I write it like this
data <- read_in_data(inputfile)
toPlot <- subset(data, ! polyaa %in% c("NQ5", "ST5") & num == 5)
gp <- plot_polydata(toPlot)

ggsave(filename = outputpdf,
       plot = gp,
       useDingbats=FALSE)