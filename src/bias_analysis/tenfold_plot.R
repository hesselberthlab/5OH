# Create boxplot of reads by ss/ds status
# Input is tab file with two columns: tenfold status (0=ss, 1=ds)
# and UMI-corrected read counts

# S3_structure is from intersection with just positive strand data

library(ggplot2)

#Read data
count_data = read.table("S3_structure.tab", header=TRUE)
dfx <- read.delim("S3_structure.tab", header=TRUE)

# Define plotting function which takes data frame and a minimum count threshold
plot_structure <- function(dframe,min_count){
  gp <- ggplot(subset(dframe, count > min_count-1), aes(factor(tenfold), count))
  gp + geom_boxplot() + 
    coord_trans(y = "log10") +
    labs(title=paste("Reads by secondary structure prediction\nMin count: ", min_count),
         x="Single stranded                 Double stranded",
         y="UMI Reads")
}

plot_structure(dfx,0)
plot_structure(dfx,50)
plot_structure(dfx,100)

