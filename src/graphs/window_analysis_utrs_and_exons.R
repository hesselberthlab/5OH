# Middle version for exon + UTR graphs.

library(dplyr)
library(ggplot2)
library(reshape2)

# import data

win10 <- read.table("ATTGGC_S4.window.full.tab",header=TRUE) %.%
  select(gene,bin,count,cat,strand) %.%
  filter(count>9)
win10$gene <- factor(win10$gene)

win40 <- read.table("ATTGGC_S4.window.full.tab",header=TRUE) %.%
  select(gene,bin,count,cat,strand) %.%
  filter(count>39)
win40$gene <- factor(win40$gene)

win50 <- read.table("ATTGGC_S4.window.full.tab",header=TRUE) %.%
  select(gene,bin,count,cat,strand) %.%
  filter(count>49)
win50$gene <- factor(win50$gene)

difdisc <- c("RPS31","TPI1","ADE8","TDH1","MDH1","RPN2")


# get window data, sum counts in same gene bin
win <- read.table("ATTGGC_S4.window.full.tab",header=TRUE) %.%
  select(gene,bin,count,cat,strand) %.%
  filter(count>9)
  #filter(gene %in% difdisc, count>24)
  #filter(gene %in% win50$gene)#,count>4)
  #filter(! gene %in% win10$gene)
win$bin <- as.integer(as.character(win$bin))
win <- win %.% 
  group_by(gene,bin) %.% 
  summarize(counts=sum(count)) %.%
  #summarize(counts=log2(sum(count))) %.%
  arrange(desc(counts)) # arrange by counts
  #arrange(desc(bin),desc(counts))  # arrange by bins and by count

win$gene_names <- factor(win$gene, as.character(win$gene)) # force R to plot in this order

# scatter plot of counts vs. bin
#biggies <- win %.% filter(counts>2000)
#gp <- ggplot(win,aes(x=bin,y=counts,colour=counts))
#gp + geom_point() + 
#  geom_text(data=biggies,aes(label=gene),size=4,hjust=1.1) +
#  theme_bw()
#scale_colour_manual(values = c("-"="red", "+"="blue"))


# Heatmap-esque graph
base_size <- 9

UTR_3 <- data.frame(x = 0.5, y =c(-Inf,Inf),UTR_3 = factor("3' UTR") )
UTR_5 <- data.frame(x = 20.5, y =c(-Inf,Inf),UTR_5 = factor("5' UTR") )
#middle <- data.frame(x = 10.5, y =c(-Inf,Inf), middle = factor("middle") )

hm <- win #%.% filter(counts>5)

to_label <- subset(hm, gene %in% c("ADE8","RPS31","TPI1","MDH1","TDH1","SPT6"))

#p <- ggplot(hm, aes(y=bin,x=gene_names))
p <- ggplot(hm, aes(y=reorder(gene_names,bin),x=bin))
gene_num <- win %.% group_by(gene) %.% summarize(n()) %.% nrow()
p + geom_tile(aes(fill = counts)) +
  scale_fill_gradient() +
  theme_bw() + 
  labs(x = "5'UTR         CDS          3'UTR",y = paste("genes (", as.character(gene_num),")",sep="")) + 
  geom_line(aes(x, y,linetype = "solid"), UTR_3, show_guide=FALSE) + 
  geom_line(aes(x, y,linetype = "solid"), UTR_5, show_guide=FALSE) +
  #geom_text(data=to_label,aes(label=gene,y=10),angle=90,size=3) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(limit = c(-1, 22)) +
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.x = element_text(angle=90,size=7),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle=90))