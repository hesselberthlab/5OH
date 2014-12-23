# Simple version for exon-only graphs.  Totally works.

library(dplyr)
library(ggplot2)
#library(reshape2)

# import data

# get window data, sum counts in same gene bin
win <- read.table("ATTGGC_S4.window.bg",header=TRUE) %.%
  select(gene,bin,count) %.%
  filter(count>9)
win$bin <- as.integer(as.character(win$bin))
win <- win %.% group_by(gene,bin) %.% summarize(counts=sum(count))

# scatter plot of counts vs. bin
biggies <- win %.% filter(counts>2000)
gp <- ggplot(win,aes(x=bin,y=counts,colour=counts))
gp + geom_point() + 
  geom_text(data=biggies,aes(label=gene),size=4,hjust=1.1) +
  theme_bw()
  #scale_colour_manual(values = c("-"="red", "+"="blue"))


# heatmap attempts...
hm <- win %.% filter(counts<3000,counts>199)

base_size <- 9

p <- ggplot(hm, aes(x=bin,y=reorder(gene,bin)))
p + geom_tile(aes(fill = counts),colour = "white") +
  scale_fill_gradient(low = "black", high = "red") +
  theme_grey(base_size = base_size) + 
  labs(x = "bins",y = "genes") + 
  #scale_y_reverse() +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none",
       axis.ticks = element_blank(),
       axis.text.x = element_text(size = base_size*0.8,hjust = 0, 
                                colour = "grey50"))
