import dplyr
import reshape2
import ggplot2

# Parse Hammarlund data and generate scatter plot by gene

gene_exp <- read.table("gene_exp_diff.tab", header=FALSE)
colnames(gene_exp) <- c("test_id", "gene_id", "gene", "locus",
                        "sample_1", "sample_2",  "status",
                        "value_1", "value_2", "log2(fold_change)",
                        "test_stat", "p_value", "q_value", "significant", "gene_name")
df_1 <- gene_exp %.% select(gene_name,sample_1,value_1) %.% unique()
df_2 <- gene_exp %.% select(gene_name,sample_2,value_2) %.% unique()
colnames(df_2) <- c("gene_name","sample_1","value_1")

plot_data <- merge(df_1,df_2,all=TRUE) %.%
  filter(value_1 > 1) %.%
  group_by(gene_name) %.%
  filter(n()>3)

plot_data$sample_1 <- factor(plot_data$sample_1,
                      levels=c("CNTRL", "CNTRL_TM", "RTCB","RTCB_TM"),
                      labels=c("Control", "Control + Tm",
                               "RtcB","RtcB + Tm"))

gp <- ggplot(plot_data,aes(x=gene_name,y=value_1,colour=sample_1))
gp + geom_point(size=3.5) +
  ylab("RPKM") +
  xlab("Gene") +
  scale_y_log10() +
  labs(colour="Treatment") +
  theme(axis.text.x = element_text(size=12,angle=45, vjust=0.75, color="black"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12))