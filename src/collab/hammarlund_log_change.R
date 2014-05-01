import dplyr
import reshape2
import ggplot2

# Parse Hammarlund RPKM expression data and generate frequency plot by gene

gene_exp <- read.table("gene_exp.diff", header=TRUE)

ctrl_ind <- gene_exp %.%
  select(test_id,sample_1,sample_2,value_1,value_2,log2.fold_change.) %.%
  filter(sample_1=="CNTRL",
         sample_2=="CNTRL_TM",
         log2.fold_change.>=1.75,
         value_1>=1,
         value_2>=1)

rtcb_ind <- gene_exp %.%
  select(test_id,sample_1,sample_2,value_1,value_2,log2.fold_change.) %.%
  filter(sample_1=="RTCB",
         sample_2=="RTCB_TM",
         #log2.fold_change.>=1.75,
         value_1>=1,
         value_2>=1)

# _02 graph 
plot_data <- merge(ctrl_ind,rtcb_ind,by="test_id") %.%
  arrange(desc(log2.fold_change..y))

# _01 graph
plot_data <- merge(ctrl_ind,rtcb_ind,by="test_id")  %.%
  filter(#value_1.x <= value_1.y*100,
         #value_1.x >= value_1.y/100,
         value_1.y <= value_1.x) %.%
  arrange(desc(log2.fold_change..x))

gp <- ggplot(plot_data,aes(x=rev(reorder(test_id, log2.fold_change..y)),y=log2.fold_change..x))
gp + geom_bar(aes(fill="Control"),stat="identity") +
  geom_bar(data=plot_data,stat="identity",aes(y=log2.fold_change..y,fill="RtcB"),alpha=0.5) +
  ylab("Log2 Fold Change, Tm vs. DMSO") +
  xlab("Gene") +
  labs(fill="Strain") +
  scale_fill_manual(values=c("lightblue", "#CC79A7"), 
                    name="Strain",
                    labels=c("Control", "RtcB")) +
  scale_y_continuous(breaks=c(-2,0,2,4)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size=10,angle=90,vjust=0.7, color="black"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))