# discrete_hits_identification.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# Identify discrete hits in 5OH-Seq data;
# ie, sites with signal at least 3-fold higher than the second highest signal
# Compare align-all with align-uniq for WT_DMSO sample.

library(dplyr)

# Generates list of most abundant gene site, the proportion of signal in that site, the next highest proportion
# and the fold increase in proportion in the highest signal vs. the second highest.
# list is sorted by fold increase, such that the most unique cleavages are at the top of the list
get_single_hits <- function(df) {
  comp.props <- data.frame()
  for (gene.name in levels(factor(df$gene))) {
    maxgene <- subset(df, gene == gene.name) %>%
      arrange(desc(proportion)) %>%
      mutate(site = paste(chrm, start, sep=":"))
    gene.site <- maxgene[1,]$site
    site.prop <- maxgene[1,]$proportion
    next.prop <- maxgene[2,]$proportion
    fold.incr <- site.prop / next.prop
    comp.props <- rbind(comp.props, data.frame(gene.name, gene.site, site.prop, next.prop, fold.incr))
  }
  comp.props <- comp.props %>%
    arrange(desc(fold.incr)) %>%
    subset(fold.incr >= 3.0)
  return(comp.props)
}

get_prop_comp <- function(sample.name, assembly, alignment) {
  fileprefix <- paste(sample.name, assembly, alignment, sep=".")
  table <- paste("proportions/", fileprefix, ".props.bed", sep="")
  df <- read.table(table, header=FALSE, col.names = c("chrm", "start", "stop", "gene", "proportion", "strand"))
  comp.props <- get_single_hits(df)
  return(comp.props)
}
  
}

# Alignment specifics
workingdir = "~/projects/5OH/results/methodspaper/"
setwd(workingdir)
samples <- c("SP8")
assembly <- "assembly-sacCer1"
alignments <- c("align-uniq", "align-all")

# Parse input file
for (sample in samples) {
  for (alignment in alignments) {
    comp.props <- get_prop_comp(sample, assembly, alignment)
    print(comp.props)
  }
}


############################################################
# Comparison of align-all and align-uniq for single sample #
############################################################

sample = "SP8"
assembly <- "assembly-sacCer1"

# Align all
alignment <- "align-all"
alignall.comp.props <- get_prop_comp(sample, assembly, alignment)
alignall.genes <- data.frame(levels(factor(alignall.comp.props$gene.name)))
colnames(alignall.genes) <- c("gene.name")

# Align uniq
alignment <- "align-uniq"
alignuniq.comp.props <- get_prop_comp(sample, assembly, alignment)
alignuniq.genes <- data.frame(levels(factor(alignuniq.comp.props$gene.name)))
colnames(alignuniq.genes) <- c("gene.name")

# Which of the align-uniq hits are also represented in the align-all set?
subset(alignuniq.genes, gene.name %in% alignall.genes$gene.name)

# Which of the align-uniq hits are not also represented in the align-all set?
subset(alignuniq.genes, ! gene.name %in% alignall.genes$gene.name)

# Which of the align-all hits are not represented in the align-uniq set?
subset(alignall.comp.props, ! gene.name %in% alignuniq.genes$gene.name)