genedf <- subset(df, gene == "YLR167W")
m1 <- glm.nb(count.level ~ 1, data=genedf)
sites <- nrow(genedf)
badasstable <- genedf %>%
  mutate(p.value = 1 - (pnbinom(count.level,m1$theta,mu=m1$fit[1]))) %>%
  mutate(BH.p.value = p.adjust(p.value, method = "BH", n=sites)) %>%
  arrange(p.value)
badasstable %>% head(10)

pvals <- badasstable$p.value
qvals <- p.adjust(pvals, method = "BH")
qvalue(pvals)