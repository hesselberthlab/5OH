library(ggplot2)
library(dplyr)

colnames <- c('pos','rel.pos','signal','library.type','feature.label',
              'sample.label')
df <- read.table('combined.tab', col.names=colnames)

gp <- ggplot(aes(x = rel.pos, y = signal, color=factor(sample.label)), data = df)

gp <- gp + geom_line(size=1) 
gp <- gp + scale_color_brewer(palette="Set1")
gp <- gp + geom_vline(xintercept=0)
gp <- gp + theme_bw()

gp <- gp + ggtitle('hClp1 KD density plots')
gp <- gp + xlab('Relative position (bp)')
gp <- gp + ylab('Signal (sum of means)')

gp <- gp + theme(legend.position = 'bottom')
gp <- gp + guides(color = guide_legend(title = NULL))

gp <- gp + facet_grid(library.type ~ feature.label, scales="free")


gp
