import dplyr
import reshape2
import ggplot2

# Parse Hammarlund data and generate scatter plot by gene

gene_exp <- read.table("gene_exp_diff.tab", header=FALSE)
colnames(gene_exp) <- c("test_id", "gene_id", "gene", "locus",
                        "sample_1", "sample_2",  "status",
                        "value_1", "value_2", "log2(fold_change)",
                        "test_stat", "p_value", "q_value", "significant")
df_1 <- gene_exp %.% select(gene_id,sample_1,value_1) %.% unique()
df_2 <- gene_exp %.% select(gene_id,sample_2,value_2) %.% unique()
colnames(df_2) <- c("gene_id","sample_1","value_1")

plot_data <- merge(df_1,df_2,all=TRUE) %.%
  filter(value_1 > 0)

gp <- ggplot(plot_data,aes(x=gene_id,y=value_1,colour=sample_1))
gp + geom_point(size=3.5) +
  ylab("RPKM") +
  xlab("Gene") +
  labs(colour="Treatment") +
  scale_y_log10()