import dplyr
import reshape2
import ggplot2

# Parse Hammarlund RPKM expression data and generate boxplots by gene

gene_exp <- read.table("gene_exp.diff", header=TRUE)
df_1 <- gene_exp %.% select(test_id,sample_1,value_1) %.% unique()
df_2 <- gene_exp %.% select(test_id,sample_2,value_2) %.% unique()
colnames(df_2) <- c("test_id","sample_1","value_1")

plot_data <- merge(df_1,df_2,all=TRUE) %.%
  filter(value_1 >= 1) %.%
  group_by(test_id) %.%
  filter(n()>3)

plot_data$sample_1 <- factor(plot_data$sample_1,
                             levels=c("CNTRL", "CNTRL_TM", "RTCB","RTCB_TM"),
                             labels=c("Control", "Control + Tm",
                                      "RtcB","RtcB + Tm"))

gp <- ggplot(plot_data,aes(x=sample_1,y=value_1,fill=sample_1))
gp + geom_boxplot(color="black") +
  ylab("RPKM") +
  xlab("Treatment") +
  scale_y_log10(breaks=c(10, 100, 1000, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12,vjust=0.75, color="black"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.position="none")

print ("Number of Genes:")
print (nrow(plot_data)/4)