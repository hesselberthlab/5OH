# histogram of count frequencies

library(dplyr)
library(ggplot2)
library(reshape2)

wt <- read.table("ATTGGC_S4.bg",header=FALSE)
tm <- read.table("CACTGT_S3.bg",header=FALSE)
names(wt) = c("chrm","start","stop","reads")
names(tm) = c("chrm","start","stop","reads")

wt <- subset(wt,reads<1500) %.% mutate(sample="WT_DMSO")
tm <- subset(tm,reads<1500) %.% mutate(sample="WT_tm")

samples <- merge(wt,tm,all=TRUE)

ggplot(samples,aes(x=sample,y=reads)) +
  geom_boxplot() + geom_jitter()

ggplot(wt, aes(x=reads)) +
  geom_histogram(binwidth=25) + 
  scale_y_log10()

trimcats <- c("gene","CDS_Verified","noncoding_exon",
              "CDS_Dubious","CDS_Uncharacterized",
              "CDS_Verified%7Csilenced_gene","pseudogene",
              "gene_cassette","Y_prime_element",
              "internal_transcribed_spacer_region",
              "external_transcribed_spacer_region",
              "nucleotide_match","non_transcribed_region",
              "transposable_element_gene","X_element_combinatorial_repeat",
              "centromere","centromere_DNA_Element_II","five_prime_UTR_intron",
              "centromere_DNA_Element_III", "X_element","CDS",
              "intron","long_terminal_repeat","LTR_retrotransposon",
              "telomere")

### Loading data; write better code for such things.
tm <- read.table("CACTGT_S3.intersect.tab",header=TRUE) %.%
  mutate(sample="WT_Tm_rep1",strain="WT",condition="Tun",SAP="-")
wt <- read.table("ATTGGC_S4.intersect.tab",header=TRUE)  %.%
  mutate(sample="WT_DMSO_rep1",strain="WT",condition="DMSO",SAP="-")
tm_sap <- read.table("AAGCTA_S7.intersect.tab",header=TRUE)  %.%
  mutate(sample="WT_Tm_SAP_rep1",strain="WT",condition="Tun",SAP="+")
wt_sap <- read.table("GTAGCC_S8.intersect.tab",header=TRUE)  %.%
  mutate(sample="WT_DMSO_SAP_rep1",strain="WT",condition="DMSO",SAP="+")

xtm <- read.table("GCCTAA_S1.intersect.tab",header=TRUE) %.%
  mutate(sample="Xrn1_Tm_rep1",strain="Xrn1",condition="Tun",SAP="-")
xwt <- read.table("TGGTCA_S2.intersect.tab",header=TRUE)  %.%
  mutate(sample="Xrn1_DMSO_rep1",strain="Xrn1",condition="DMSO",SAP="-")
xtm_sap <- read.table("TCAAGT_S5.intersect.tab",header=TRUE)  %.%
  mutate(sample="Xrn1_Tm_SAP_rep1",strain="Xrn1",condition="Tun",SAP="+")
xwt_sap <- read.table("CTGATC_S6.intersect.tab",header=TRUE)  %.%
  mutate(sample="Xrn1_DMSO_SAP_rep1",strain="Xrn1",condition="DMSO",SAP="+")

wtdmso <- merge(wt,tm,all=TRUE)
wtsap <- merge(wt_sap,tm_sap,all=TRUE)
wt_all <- merge(wtdmso,wtsap,all=TRUE)

xdmso <- merge(xwt,xtm,all=TRUE)
xsap <- merge(xwt_sap,xtm_sap,all=TRUE)
x_all <- merge(xdmso,xsap,all=TRUE)

all_samples <- merge(wt_all,x_all,all=TRUE)


# Graph RNAs by chromosome; ~not particularly useful
ggplot(merge(wt,tm,all=TRUE),aes(x=chr, y=count,fill=cat,colour=cat)) +
  geom_boxplot() + geom_jitter() + facet_grid(~sample) +
  labs(x="Chromosomes",y="Counts Per Millions (CPM)",title="Types of RNA by Chromosome") +
  scale_colour_discrete(guide=FALSE) +
  scale_x_discrete(breaks=NULL) +
  guides(fill=guide_legend(title="Types Of RNA"))

# Graph RNAs by conditions, color by type
# change subset argument to yield different plots.
ggplot(subset(wt_all,! cat %in% trimcats),aes(x=condition, y=count,fill=cat,colour=cat)) +
  geom_boxplot(outlier.size=0) + geom_jitter(size=2) + facet_grid(~SAP) +
  labs(x="Samples",y="Counts Per Millions (CPM)",title="Sites Identified by 5OH Assay") +
  scale_colour_discrete(guide=FALSE) +
  guides(fill=guide_legend(title="Types Of RNA"))

data = subset(wt_all,count > 10 & cat == "intron" & gene != "EFB1")

ggplot(data,aes(x=SAP, y=count,fill=cat,colour=cat)) +
  geom_boxplot(outlier.size=0) + geom_jitter(size=2) + facet_grid(~gene) +
  labs(x="Samples",y="Counts Per Millions (CPM)",title="Types of RNA by Sample") +
  scale_colour_discrete(guide=FALSE) +
  guides(fill=guide_legend(title="Types Of RNA")) + 
  geom_text(data=data, mapping=aes(label=gene),
            size=3.5,hjust=1.1, vjust=0)


## analyze top genes...


##
wt2 <- wt %.%
  mutate(wt.ct=count) %.%
  filter(cat=="mRNA") %.%
  select(gene,site,wt.ct)
tm2 <- tm %.%
  mutate(tm.ct=count) %.%
  filter(cat=="mRNA") %.%
  select(gene,site,tm.ct)

wt_trmt <- inner_join(wt2,tm2)

wt_labels <- subset(wt_trmt, wt.ct >= tm.ct * 2)
tm_labels <- subset(wt_trmt, tm.ct >= wt.ct * 2)

ggplot(wt_trmt, aes(x=wt.ct,y=tm.ct)) +
  geom_point() +
  xlim(0,500) + ylim(0,500) +
  geom_smooth(method = "lm", color="black", formula = y ~ x) +
  geom_point(data=wt_labels,color="red") +
  geom_point(data=tm_labels,color="darkgreen") +
  geom_text(data=filter(wt_labels,wt.ct > 150),mapping=aes(label=gene),size=4) + 
  geom_text(data=filter(tm_labels,tm.ct > 300),mapping=aes(label=gene),size=4)
  