# Work in progress to identify cleavage changes in experimental conditions
# need to turn everything into functions...

library(dplyr)
library(ggplot2)
library(reshape2)

tm <- read.table("CACTGT_S3.intersect.tab",header=TRUE) %.%
  mutate(sample="WT_Tm_rep1",strain="WT",condition="Tun",SAP="-")
wt <- read.table("ATTGGC_S4.intersect.tab",header=TRUE)  %.%
  mutate(sample="WT_DMSO_rep1",strain="WT",condition="DMSO",SAP="-")
xwt <- read.table("TGGTCA_S2.intersect.tab",header=TRUE)  %.%
  mutate(sample="Xrn1_DMSO_rep1",strain="Xrn1",condition="DMSO",SAP="-")

### Original raw code
wt_mRNA <- wt %.%
  mutate(wt.ct=count) %.%
  filter(cat=="mRNA") %.%
  select(gene,site,wt.ct)
tm_mRNA <- tm %.%
  mutate(tm.ct=count) %.%
  filter(cat=="mRNA") %.%
  select(gene,site,tm.ct)
wt_trmt <- inner_join(wt_mRNA,tm_mRNA)

wt_trmt2 <- inner_join(get_mRNA(wt),get_mRNA(tm),by=c("gene","site"))
wt_change <- subset(wt_trmt, wt.ct >= tm.ct * 2)
tm_change <- subset(wt_trmt, tm.ct >= wt.ct * 2)

ggplot(wt_trmt, aes(x=wt.ct,y=tm.ct)) +
  geom_point() +
  xlim(0,500) + ylim(0,500) +
  geom_smooth(method = "lm", color="black", formula = y ~ x) +
  geom_point(data=wt_change,color="red") +
  geom_point(data=tm_change,color="darkgreen") +
  geom_text(data=filter(wt_change,wt.ct > 150),mapping=aes(label=gene),size=4) + 
  geom_text(data=filter(tm_change,tm.ct > 300),mapping=aes(label=gene),size=4)

## Sort sites by fold change; num rows = num sites
tm_change %.% filter(tm.ct > 39) %.% group_by(gene) %.% arrange(desc(tm.ct/wt.ct))

## Identify number of genes with 2-fold change; num rows = num genes
tm_change %.% filter(tm.ct > 39) %.% group_by(gene) %.% summarize(n())

tm_change %.% filter(tm.ct > 39) %.% group_by(gene) %.% arrange(desc(tm.ct))


### Functionalization of code
get_mRNA <- function(df) {
  df <- df %.%
    filter(cat=="mRNA") %.%
    select(gene,site,count)
  
  return(df)
}

plot_site_counts <- function(df,x_changes,y_changes,dfx.label,dfy.label) {
  gp <- ggplot(df, aes(x=count.x,y=count.y)) +
    geom_point() +
    xlim(0,500) + ylim(0,500) +
    geom_smooth(method = "lm", color="black", formula = y ~ x) +
    geom_point(data=x_changes,color="red") +
    geom_point(data=y_changes,color="darkgreen") +
    xlab(dfx.label) + ylab(dfy.label) + 
    geom_text(data=filter(x_changes,count.x > 150),mapping=aes(label=gene),size=4) + 
    geom_text(data=filter(y_changes,count.y >= count.x*10,count.y>100),
              mapping=aes(label=gene),size=4,
              vjust=1)
  
  gp
}

stats_for_site_comparison <- function(x,y,dfx.label,dfy.label) {
  x_sites <- x %.% 
    filter(count.x > 39) %.% 
    group_by(gene) %.% 
    arrange(desc(count.x/count.y))
  x_genes <- x_sites %.%
    summarize(n())
  
  y_sites <- y %.% 
    filter(count.y > 39) %.% 
    group_by(gene) %.% 
    arrange(desc(count.y/count.x))
  y_genes <- y_sites %.%
    summarize(n())

  print (paste("Num sites increased in ",dfx.label,": ",nrow(x_sites),sep=""))
  print (paste("Num genes increased in ",dfx.label,": ",nrow(x_genes),sep=""))
  
  print (paste("Num genes increased in ",dfy.label,": ",nrow(y_genes),sep=""))
  print (paste("Num sites increased in ",dfy.label,": ",nrow(y_sites),sep=""))
  
  return (list(x=x_sites,y=y_sites))
}

site_plot <- function(dfx,dfy,dfx.label,dfy.label,fold) {
  dfx_mRNA <- get_mRNA(dfx)
  dfy_mRNA <- get_mRNA(dfy)
  
  dfs_join <- inner_join(dfx_mRNA,dfy_mRNA, by=c("gene","site"))
  x_change <- subset(dfs_join, count.x >= count.y * fold)
  y_change <- subset(dfs_join, count.y >= count.x * fold)
  
  plot_site_counts(dfs_join,x_change,y_change,dfx.label,dfy.label)
}

site_stats <- function(dfx,dfy,dfx.label,dfy.label,fold) {
  dfx_mRNA <- get_mRNA(dfx)
  dfy_mRNA <- get_mRNA(dfy)
  
  dfs_join <- inner_join(dfx_mRNA,dfy_mRNA, by=c("gene","site"))
  x_change <- subset(dfs_join, count.x >= count.y * fold)
  y_change <- subset(dfs_join, count.y >= count.x * fold)

  stats_for_site_comparison(x_change,y_change,dfx.label,dfy.label)
}

site_plot(wt,wt,"DMSO","DMSO",2)
site_plot(wt,tm,"DMSO","Tm",3)
site_stats(wt,tm,"DMSO","Tm",10)

site_plot(wt,xwt,"DMSO","del_Xrn1",2)
site_stats(wt,xwt,"DMSO","del_Xrn1",5)







