library(ggplot2)
library(reshape2)
library(stringr)
library(Cairo)

# some initial sanity check graphs
setwd("~/projects/5OH/data/GeoSubmission/Rep1/etc")
plot_singlepolyaa <- function(filename,main_title) {
  data <- read.table(filename, header=FALSE)
  plot(apply(data,2,mean), main=main_title, ylab="avg CPMs", xlab="Index")
}

plot_singlepolyaa("polyEmatrix.tab","EEEEE")
plot_singlepolyaa("polyRmatrix.tab","RRRRR")
plot_singlepolyaa("polyAmatrix.tab","AAAAA")
plot_singlepolyaa("polyDmatrix.tab","DDDDD")
plot_singlepolyaa("polyKmatrix.tab","KKKKK")
plot_singlepolyaa("polyNmatrix.tab","NNNNN")


# Full data graphs and subset graphs
read_in_data <- function(filename) {
  data <- read.csv(filename)
  data <- melt(data, id.vars = c("bin"),
               variable.name = "polyaa")
  data$value <- as.numeric(data$value)
  
  parsed_polyaa = as.data.frame(str_match(data$polyaa, "^([A-Z]|[A-Z][A-Z]|[A-Z][A-Z][A-Z])([1-9]|[1-9][0-9])$")[,-1])
  colnames(parsed_polyaa) <- c("aa", "num")
  data <- cbind(data, parsed_polyaa)
  
  return(data)
}

plot_polydata_subsets <- function(df,polyaa_list) {
  ggplot(subset(data, polyaa %in% polyaa_list), aes(x=index, y=value, color=polyaa)) +
    geom_point() +
    theme_bw() +
    ylab("Proportion x 10,000") +
    scale_x_continuous(breaks=c(0,25,50,75,100), labels=c("-50", "-25", "0", "25", "50")) +
    xlab("Distance from poly stretch (bp)")
}

plot_polydata <- function(df) {
  ggplot(data, aes(x=bin, y=log10(value), color=polyaa)) +
    geom_point() +
    theme_bw() +
    ylab("Proportion x 10,000") +
    scale_x_continuous(breaks=c(0,25,50,75,100), labels=c("-50", "-25", "0", "25", "50")) +
    xlab("Distance from poly stretch (bp)") +
    facet_wrap(~aa) +
    theme(legend.position="none")
}

# WT DMSO data
setwd("~/projects/5OH/data/GeoSubmission/Rep1/etc")
data <- read_in_data("WT.DMSO.Rep1.poly-all.avgs.csv")
plot_polydata(data)

# Some interesting subsets
polyfive <- c("A5", "D5", "E5", "H5", "I5", "K5", "L5", "N5", "P5", "Q5", "R5", "S5", "T5")
polyE <- c("E5", "E10", "E15")
polyK <- c("K5", "K10", "K15")
polyregex <- c("DE.5", "DE.10", "DE.15", "RK.5", "RK.10", "RL.15")
polybasic <- c("R5", "K5", "RK5")
polyacidic <- c("D5", "E5", "DE5")
polyeracidic <- c("DE5", "DE10", "DE15")

# where the magic happens
plot_polydata_subsets(data,polyfive)
plot_polydata_subsets(data,polybasic)
plot_polydata_subsets(data,polyacidic)
plot_polydata_subsets(data,polyeracidic)

# Ribozero data
setwd("~/projects/5OH/results/20140529_ctrls")
data <- read_in_data("ribozero.poly-all.avgs.csv")
plot_polydata(data)

# Roy Parker data
setwd("~/projects/collab/parker")
data <- read_in_data("delXrn1Dcp2.poly-all.avgs.csv")
plot_polydata(data)
plot_polydata_subsets(data,polybasic)
plot_polydata_subsets(data,polyacidic)


################################
# Merged and colored by sample #
################################

plot_polydata_merge <- function(df) {
  ggplot(data, aes(x=index, y=value, color=sample)) +
    geom_point() +
    theme_bw() +
    ylab("Average CPM") +
    scale_x_continuous(breaks=c(0,25,50,75,100), labels=c("-50", "-25", "0", "25", "50")) +
    xlab("Distance from poly stretch (bp)") +
    facet_wrap(~aa) +
    theme(legend.position="none")
}

# Roy Parker data
setwd("~/projects/collab/parker")
Pdata <- read_in_data("delXrn1Dcp2.poly-all.avgs.csv") %.% mutate(sample = "5P")

# WT DMSO data
setwd("~/projects/5OH/data/GeoSubmission/Rep1/etc")
OHdata <- read_in_data("WT.DMSO.Rep1.poly-all.avgs.csv") %.% mutate(sample = "5OH")

plot_polydata_merge(rbind(Pdata, OHdata))


###################################
# Plotting by coverage proportion #
###################################
setwd("~/projects/5OH/data/GeoSubmission/Rep1/etc")
data <- read_in_data("WT.DMSO.Rep1.props.poly-all.avgs.csv")

gp <- plot_polydata(data)
ggsave(filename = "WT.DMSO.Rep1.props-5bins-avg.poly-all.pdf", 
       plot = gp,
       device = CairoPDF)

gp <- plot_polydata_subsets(data,polybasic)
ggsave(filename = "WT.DMSO.Rep1.props.polybasic.pdf", 
       plot = gp,
       device = CairoPDF)

gp <- plot_polydata_subsets(data,polyacidic)
ggsave(filename = "WT.DMSO.Rep1.props.polyacidic.pdf", 
       plot = gp,
       device = CairoPDF)
