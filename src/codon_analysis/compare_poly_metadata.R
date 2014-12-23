# Compare poly_metadata between samples

library(ggplot2)
library(dplyr)
library(reshape2)

read_in_data <- function(filename) {
  data <- read.csv(filename)
  data <- melt(data, id.vars = c("index"),
               variable.name = "polyaa")
  data$value <- as.numeric(data$value)
  
  parsed_polyaa = as.data.frame(str_match(data$polyaa, "^([A-Z]|[A-Z][A-Z]|[A-Z][A-Z][A-Z])([1-9]|[1-9][0-9])$")[,-1])
  colnames(parsed_polyaa) <- c("aa", "num")
  data <- cbind(data, parsed_polyaa)
  
  return(data)
}

get_ratio <- function(df,num_bp) {
  data25 <- subset(df, index==50 - num_bp)
  colnames(data25)[3] <- "value.25"
  data75 <- subset(df, index==50 + num_bp)
  colnames(data75)[3] <- "value.75"
  merged_df <- cbind(data25, data75[3]) %.%
    mutate(ratio = value.25 / value.75)
  
  return(subset(merged_df, value.25 > 0 & value.75 >0))
}

# WT DMSO data
setwd("~/projects/5OH/data/GeoSubmission/Rep1/etc")
OHdata <- read_in_data("WT.DMSO.Rep1.poly-all.avgs.csv") %.% mutate(sample = "5OH")

# Roy Parker data
setwd("~/projects/collab/parker")
Pdata <- read_in_data("delXrn1Dcp2.poly-all.avgs.csv") %.% mutate(sample = "5P")

# Merge data, calculate ratios upstream and downstream of polyaa motif start, and plot
alldata <- rbind(OHdata, Pdata)
ratios <- get_ratio(alldata, 25)
aa_motifs <- subset(ratios, polyaa %in% c("K5", "R5", "RK5", "E5", "D5", "DE5"))
ggplot(aa_motifs, aes(x=polyaa, y=ratio, fill=sample)) +
  geom_bar(stat="identity", position="dodge") +
  ylab("Ratio: (Signal at -25) / (Signal at 25)") +
  xlab("PolyAA motif")

