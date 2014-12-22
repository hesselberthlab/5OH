library(dplyr)
library(reshape2)
library(ggplot2)

# Parse Hammarlund RPKM expression data and generate frequency plot by gene

gene_exp <- read.table("gene_exp.diff", header=TRUE)
upr_genes <- read.table("upr_targets.tab", header=TRUE)
colnames(upr_genes) <- c("test_id","gene","gene_id")
upr_ids <- upr_genes$test_id



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
         #log2.fold_change.>=1.75, # for _03 graph, remove comment
         value_1>=1,
         value_2>=1)

# _02 graph 
plot_data <- merge(ctrl_ind,rtcb_ind,by="test_id") %.%
  arrange(desc(log2.fold_change..x))

plot_data$test_id = factor(plot_data$test_id, levels=reorder(plot_data$test_id, plot_data$log2.fold_change..x))

gene_ids <- plot_data$test_id
names_to_plot <- plot_data %.% filter(test_id %in% upr_ids)
names_to_plot <- merge(names_to_plot,upr_genes,by.x="test_id")
names_to_plot <- names_to_plot %.% select(test_id)

# _01 graph
#plot_data <- merge(ctrl_ind,rtcb_ind,by="test_id")  %.%
#  filter(#value_1.x <= value_1.y*100,
#         #value_1.x >= value_1.y/100,
#         value_1.y <= value_2.x)

graph_control_desc(plot_data)
#graph_rtcb_desc(plot_data)

graph_control_desc <- function(df) {
  gp <- ggplot(df,aes(x=test_id,width=0.85))
  gp + 
    geom_bar(aes(y=log2.fold_change..x,fill="Control"),stat="identity") +
    geom_bar(aes(y=log2.fold_change..y,fill="RtcB"),stat="identity",alpha=0.5) +
    geom_bar(data=names_to_plot,stat="identity",aes(x=test_id,
                                                    y=log2.fold_change..y,
                                                    fill="Known targets"),alpha=0.5) +
    geom_text(data=names_to_plot,aes(x=test_id,label=gene_id,y=log2.fold_change..y),hjust=0,angle=90) +
    ylab("Log2 Fold Change, Tm vs. DMSO") +
    xlab("Gene") +
    labs(fill="") +
    scale_fill_manual(values=c("lightblue", "black","#CC79A7"), 
                      name="",
                      labels=c("Control", "Known XBP1 Targets","RtcB")) +
    scale_y_continuous(breaks=c(-2,0,2,4)) + 
    scale_x_discrete(breaks = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10,angle=90,vjust=0.7, color="black"),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"))
}


graph_rtcb_desc <- function(data) {
  plot_data <- data %.%
    arrange(desc(log2.fold_change..y))
  gp <- ggplot(plot_data,aes(x=rev(reorder(test_id, log2.fold_change..y)),
                             y=log2.fold_change..x,
                             width=0.85))
  gp + geom_bar(aes(fill="Control"),stat="identity") +
    geom_bar(data=plot_data,stat="identity",aes(y=log2.fold_change..y,fill="RtcB"),alpha=0.5) +
    ylab("Log2 Fold Change, Tm vs. DMSO") +
    xlab("Gene") +
    labs(fill="Strain") +
    scale_fill_manual(values=c("lightblue", "#CC79A7"), 
                      name="Strain",
                      labels=c("Control", "RtcB")) +
    scale_y_continuous(breaks=c(-4,-2,0,2,4)) + 
    scale_x_discrete(breaks = "none") +
    theme_bw() + 
    theme(#axis.text.x = element_text(size=10,angle=90,vjust=0.7, color="black"),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"))
}


