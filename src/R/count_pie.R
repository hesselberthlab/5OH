import dplyr
import reshape2
import ggplot2

# Pie chart to display number of sequencing counts
# Input unnormalized bg file
count_pie <- function(bg_filename, sample_name) {
  bg <- read.table(bg_filename,header=FALSE)
  colnames(bg) <- c("chrm","start","stop","count")
  
  one <- nrow(filter(bg,count==1))
  two_to_ten <- nrow(filter(bg,count>1,count<11))
  eleven_to_hundred <- nrow(filter(bg,count>10,count<101))
  hundredone_to_thousand <- nrow(filter(bg,count>100,count<1001))
  thousandone_plus <- nrow(filter(bg,count>1000))
  
  # If one is included
  counts <- c(one, two_to_ten, eleven_to_hundred, hundredone_to_thousand,thousandone_plus)
  ranges <- c("1","2-10","11-100","101-1000","1001+")
  sample <- c(1,1,1,1,1)
  
  # If one is *not* included...
  #counts <- c(two_to_ten, eleven_to_hundred, hundredone_to_thousand,thousandone_plus)
  #ranges <- c("2-10","11-100","101-1000","1001+")
  #sample <- c(1,1,1,1)
  
  # If <10 is *not* included...
  #counts <- c(eleven_to_hundred, hundredone_to_thousand,thousandone_plus)
  #ranges <- c("11-100","101-1000","1001+")
  #sample <- c(1,1,1)
  
  count_data <- data.frame(counts,ranges,sample)
  colnames(count_data) <- c("counts","ranges","sample")
  
  plot.title="Read Counts per 5OH Site"
  
  gp <- ggplot(data=count_data, aes(x=sample,y=counts,fill=factor(ranges)))
  gp + geom_bar(width=1,stat="identity",position="stack") +
    coord_polar(theta="y") +
    xlab('') +
    ylab('') +
    labs(fill='') +
    #labs(fill='Reads per 5OH Site') +
    ggtitle(bquote(atop(bold(.(plot.title)),
                        atop(italic(.(sample_name)), ""))))
}

count_pie("ATTGGC_S4.bg","WT DMSO Rep1")

