# window_analysis_utrs_and_exons.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# 5OH pipeline - Analysis of mRNA by 20-bin CDS and 2-bin UTR
# Plot all genes containing a 5OH site above threshold as 20-bin CDS and 2-bin UTR
# TODO: Prettify.  Create manual scale of counts.  Cull results that don't have 20-bin and 2-bin.
# Allow setting of upper threshold?

library(dplyr)
library(ggplot2)
library(reshape2)
library(Cairo)

# Parse command line arguments; end process if there is an error
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("usage: Rscript window_analysis_utrs_and_exons.R sample.windows.tab sample.name min.count output.dir")
}

sample.table = args[1]
sample.name = args[2]
min.count = as.integer(args[3])
output.dir = args[4]

# None highthruput manual entry
output.dir = "~/projects/5OH/results/methodspaper/rgraphs"
setwd(output.dir)
sample.table = "SP8.assembly-sacCer1.align-all.windows.tab"
sample.name = "WT.DMSO"
min.count = 5

# import data

# Filter data table above minimum count threshold
windows.above.threshold <- read.table(sample.table,header=TRUE) %>%
  select(gene,bin,count,cat,strand) %>%
  filter(count>=min.count)
windows.above.threshold$gene <- factor(windows.above.threshold$gene) # Refactor

# Sum counts in same gene bin
windows.above.threshold$bin <- as.integer(as.character(windows.above.threshold$bin))
windows.above.threshold <- windows.above.threshold %>% 
  group_by(gene,bin) %>% 
  summarize(counts=sum(count)) %>%
  arrange(desc(counts)) # arrange by counts
  #arrange(desc(bin),desc(counts))  # arrange by bins and by count

windows.above.threshold$gene_names <- factor(windows.above.threshold$gene,
                                             as.character(windows.above.threshold$gene)) # force R to plot in this order

# Generate graph in heatmap style

base_size <- 9
UTR_3 <- data.frame(x = 0.5, y =c(-Inf,Inf),UTR_3 = factor("3' UTR") )
UTR_5 <- data.frame(x = 20.5, y =c(-Inf,Inf),UTR_5 = factor("5' UTR") )
hm <- windows.above.threshold

p <- ggplot(hm, aes(y=reorder(gene_names,bin),x=bin))
gene_num <- windows.above.threshold %>% group_by(gene) %>% summarize(n()) %>% nrow()
p <- p + geom_tile(aes(fill = counts)) +
  scale_fill_gradient() +
  theme_bw() + 
  labs(x = "5'UTR         CDS          3'UTR",y = paste("genes (", as.character(gene_num),")",sep="")) + 
  geom_line(aes(x, y,linetype = "solid"), UTR_3, show_guide=FALSE) + 
  geom_line(aes(x, y,linetype = "solid"), UTR_5, show_guide=FALSE) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(limit = c(-1, 22)) +
  theme_bw() + 
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle=90))
p

long <- dcast(hm, gene ~ bin, value.var= "counts")

# Plot data and save file
pdf.filename <- paste(output.dir, '/', '5OH.windowanalysis.',
                      min.count, "-mincount.",
                      sample.name, '.pdf', sep='')

web.dir = "~/public_html/projects/5OH/results/methodspaper/figures"
web.filename <- paste(web.dir, '/', '5OH.windowanalysis.',
                      min.count, "-mincount.",
                      sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, 
       plot = p,
       #device = CairoPDF,
       height = 7,
       width = 4,
       useDingbats=FALSE)

ggsave(filename = web.filename, 
       plot = p,
       #device = CairoPDF,
       height = 7,
       width = 4,
       useDingbats=FALSE)