# fit_to_negative_binomial.R
#
# __author__ = 'Sally Peach'
# __contact__ = 'sallypeach@gmail.com'
#
# Fitting 5OH-Seq data to Negative Binomial Distribution
# A work in progress...

library(dplyr)
library(MASS)

results.dir <- "~/projects/5OH/results/methodspaper/proportions"
setwd(results.dir)
df <- read.table("SP8.assembly-sacCer1.align-uniq.genemap.tab.nb")

get_real_sites <- function(df, genename) {
  dftoanalyze <- subset(df, gene == genename)
  m1 <- glm.nb(count.level ~ 1, data=dftoanalyze)
  sites <- length(m1$fit) 
  
  realsitetable <- dftoanalyze %>%
    mutate(p.value = 1 - (pnbinom(count.level,m1$theta,mu=m1$fit[1]))) %>%
    mutate(BH.p.value = p.adjust(p.value, method = "BH", n=sites)) %>%
    subset(BH.p.value < FDR)
  
  return(realsitetable)
}

get_p_values <- function(df, genename) {
  dftoanalyze <- subset(df, gene == genename)
  m1 <- glm.nb(count.level ~ 1, data=dftoanalyze)
  sites <- nrow(dftoanalyze)
  
  realsitetable <- dftoanalyze %>%
    mutate(p.value = 1 - (pnbinom(count.level,m1$theta,mu=m1$fit[1]))) %>%
    mutate(BH.p.value = p.adjust(p.value, method = "BH", n=sites))
  
  return(realsitetable)
}

genes <- c("YKL085W", "YDL184C", "YIL148W")
FDR <- 0.05

corrected.sites <- data.frame()

for (genename in genes) {
  print(genename)
  realsitetable <- get_real_sites(df, genename)
  corrected.sites <- rbind(corrected.sites, realsitetable)
}

get_p_values(df, "YKL085W")

genename <- "YKL085W"
genedf <- subset(df, gene == genename)
m1 <- glm.nb(count.level ~ 1, data=dftoanalyze)
sites <- nrow(dftoanalyze)
