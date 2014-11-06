setwd("~/projects/5OH/data/GeoSubmission/Rep1/etc")
pdf.filename <- "polyaa_metadata.pdf"

data <- read.table("polyEmatrix.tab", header=FALSE)
plot(apply(data,2,mean), main="EEEEE", ylab="avg CPMs", xlab="Index")

# need better matrix skillz