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
         value_1>=1,
         value_2>=1)

p_data <- merge(ctrl_ind,rtcb_ind,by="test_id")


plot_data <- merge(ctrl_ind,rtcb_ind,by="test_id")  %.%
  filter(value_1.x <= value_1.y*10,
         value_1.x >= value_1.y/10,
         value_1.y < value_1.x) %.%
  arrange(desc(log2.fold_change..x))

gp <- ggplot(plot_data,aes(x=rev(reorder(test_id, log2.fold_change..x)),y=log2.fold_change..x))
gp + geom_bar(fill="#CC79A7",stat="identity") +
  geom_bar(data=plot_data,stat="identity",aes(y=log2.fold_change..y),fill="blue",alpha=0.4) +
  ylab("Log2 Fold Change") +
  xlab("Gene") +
  labs(fill="Strain") +
  scale_fill_manual(name="Strain",values=c("Control"="#CC79A7","RtcB"="blue"),guide="none") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=10,angle=90,vjust=0.7, color="black"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))