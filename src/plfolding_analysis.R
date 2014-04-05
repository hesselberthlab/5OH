# Create graph of read intensity vs. base pair probability
# Expectation: inverse correlation
# Actual Result: tbd...  issue comes from computing these probabilties, b/c
# inherently biased towards higher bp prob...

library(ggplot2)
library(dplyr)

pl_data = read.delim("pl_intersect.tab")
df <- subset(pl_data, count>49 & count <500)
# Subset outliers for labeling
df_o <- subset(df,count>200 & count<500)

gp <- ggplot(df, aes(x=prob,y=count))
gp +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1)) +
  labs(x="Probability of base-pairing",y="UMI-corrected Read Counts") +
  geom_text(data=df_o, mapping=aes(x=prob, y=count, label=gene),
            size=3.5,hjust=-0.1, vjust=0) +
  ggtitle(expression(atop("Read Count vs. paired probability",
                          atop(italic("WT DMSO: 20140107 S4 pos"), ""))))
ggsave("pl_intersect.png")
