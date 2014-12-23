library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(Cairo)

# some initial sanity check graphs
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

# Refactoring is the best
refactor <- function(x) {
  x <- factor(x, levels=levels(x)[levels(x) %in% x] )
  return(x)
}


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

read_in_data_mm <- function(filename) {
  data <- read.csv(filename)
  data <- melt(data, id.vars = c("index"),
               variable.name = "polyaa")
  data$value <- as.numeric(data$value)
  
  # parsing motif output, which was written as
  # AAs;ML:n;MM:n
  # ie, aminoacids of motif; motif length; mismatch number
  letters = "([A-Z]|[A-Z][A-Z]|[A-Z][A-Z][A-Z]|[A-Z][A-Z][A-Z][A-Z])"
  numbers = "([0-9]|[0-9][0-9])"
  search = paste("^", letters, "(.ML.)", numbers, "(.MM.)", numbers, "$", sep="")
  parsed_polyaa = as.data.frame(str_match(data$polyaa,search)[,-1])
  colnames(parsed_polyaa) <- c("aa", "mltxt", "motif.len", "mmtxt", "mismatch")
  parsed_polyaa <- parsed_polyaa %>%
    select(aa, motif.len, mismatch)
  data <- cbind(data, parsed_polyaa)
  
  return(data)
}

plot_polydata_subsets <- function(df,polyaa_list) {
  ggplot(subset(df, polyaa %in% polyaa_list), aes(x=index, y=value, color=polyaa)) +
    geom_point() +
    theme_bw() +
    ylab("Average CPM") +
    scale_x_continuous(breaks=c(0,25,50,75,100), labels=c("-50", "-25", "0", "25", "50")) +
    xlab("Distance from poly stretch (bp)")
}

plot_polydata <- function(df) {
  ggplot(df, aes(x=index, y=value, color=num)) +
    geom_point() +
    theme_bw() +
    ylab("Average CPM") +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75")) +
    xlab("Distance from poly stretch (bp)") +
    facet_wrap(~aa)
    #theme(legend.position="none")
}

plot_polydata_mm <- function(df) {
  ggplot(df, aes(x=index, y=value, color=mismatch)) +
    geom_point() +
    theme_bw() +
    ylab("Average CPM") +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75")) +
    xlab("Distance from poly stretch (bp)") +
    facet_wrap(~aa)
  #theme(legend.position="none")
}

plot_polydata_heatmap <- function(df) {
  df <- subset(df, num == 5) %>%
    arrange(desc(value))
  df$polyaa_sort <- factor(df$polyaa, c("RK5", "R5", "K5", "DE5", "D5", "E5", "ST5", "T5", "S5", "NQ5", "Q5", "N5", "P5", "L5", "I5", "H5", "A5"))
  df <- subset(df, polyaa_sort != "NA")
  ggplot(df, aes(x=index, y=polyaa_sort, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low="black", mid = "black", high="red", space = "rgb", name="Average\nCPMs") +
    ylab("Poly Amino Acid Motif") +
    xlab("Distance From Start of Motif (bp)") +
    theme_bw() +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75"))
}

plot_polydata_heatmap_RKED <- function(df) {
  df <- subset(df, num == 5) %>%
    arrange(desc(value))
  df$polyaa_sort <- factor(df$polyaa, c("RKED5", "RKD5","RKE5","RE5", "KE5", "RK5", "R5", "K5", "DE5", "D5", "E5", "ST5", "T5", "S5", "NQ5", "Q5", "N5", "P5", "L5", "I5", "H5", "A5"))
  df <- subset(df, polyaa_sort != "NA")
  ggplot(df, aes(x=index, y=polyaa_sort, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low="black", mid = "black", high="red", space = "rgb", name="Average\nCoverage") +
    ylab("Poly Amino Acid Motif\n(Number of motifs)") +
    xlab("Distance From Start of Motif (bp)") +
    theme_bw() +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75")) +
    scale_y_discrete(breaks=motif.labels$motif, labels=motif.labels$label)
}

plot_polydata_heatmap_nosort <- function(df) {
  #df <- subset(df, num == 5) %>%
  #  arrange(desc(value))
  #df$polyaa_sort <- factor(df$polyaa, c("KE5", "RK5", "R5", "K5", "DE5", "D5", "E5", "ST5", "T5", "S5", "NQ5", "Q5", "N5", "P5", "L5", "I5", "H5", "A5"))
  #df <- subset(df, polyaa_sort != "NA")
  ggplot(df, aes(x=index, y=polyaa, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low="black", mid = "black", high="red", space = "rgb", name="Average\nCoverage") +
    ylab("Poly Amino Acid Motif\n(Number of motifs)") +
    xlab("Distance From Start of Motif (bp)") +
    theme_bw() +
    scale_x_continuous(breaks=c(0,25,50,75,100,125,150), labels=c("-75", "-50", "-25", "0", "25", "50", "75")) +
    scale_y_discrete(breaks=motif.labels$motif, labels=motif.labels$label)
}

# WT Human data
setwd("~/projects/5OH/results/methodspaper/human")
filename <- "SP30.assembly-canonicalhg19.align-uniq.polyall.avgs.csv"
motif.labels <- read.table(paste(filename, "motifs", sep="."), header=TRUE, col.names = c("motif", "num.motifs")) %>%
  mutate(label = paste(motif, "\n", "(", num.motifs, ")", sep=""))
data <- read_in_data(filename)
plot_polydata_heatmap_RKED(data)

# WT DMSO data
setwd("~/projects/5OH/results/methodspaper")
filename <- "SP8.assembly-sacCer1.align-all.polyKRmm.avgs.csv"
data <- read_in_data_mm(filename)
plot_polydata_mm(data)

filename <- "SP8.assembly-sacCer1.align-uniq.polyall-RKED_mm2-10.avgs.csv"
data <- read_in_data_mm(filename)
plot_polydata_mm(data)

filename <- "SP8.assembly-sacCer1.align-uniq.polyall-RKED_mm3-10.avgs.csv"
data <- read_in_data_mm(filename)
motif.labels <- read.table(paste(filename, "motifs", sep="."), header=TRUE, col.names = c("motif", "num.motifs")) %>%
  mutate(label = paste(motif, "\n", "(", num.motifs, ")", sep=""))
plot_polydata_mm(data)

filename <- "SP8.assembly-sacCer1.align-uniq.polyall-RKED.avgs.csv"
motif.labels <- read.table(paste(filename, "motifs", sep="."), header=TRUE, col.names = c("motif", "num.motifs")) %>%
  mutate(label = paste(motif, "\n", "(", num.motifs, ")", sep=""))
data <- read_in_data(filename)
plot_polydata_heatmap(data)
plot_polydata_heatmap_RKED(data)
ssdata <- subset(data, aa %in% electrostatic & num == 5)
plot_polydata_heatmap_nosort(ssdata)

electrostatic <- c("R", "K", "RK", "RKE", "RKD", "RKED", "D", "E", "DE")

num.motifs <- read.table(paste(filename, "motifs", sep="."), header=TRUE) %>%
  arrange(mismatch)
num.motifs

gp <- plot_polydata_heatmap(data)

ggsave(filename = "SP8.assembly-sacCer1.align-uniq.poly-all.PolyAA-plot-motif_notworkingbooo.pdf",
       plot = gp,
       width = 8,
       height = 5)

CairoPDF(file = "WT.DMSO.Rep1.PolyAA-plot-motif.pdf",
         width = 6,
         height = 3)



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
