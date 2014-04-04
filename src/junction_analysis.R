# Creates simple R plots from junction_analysis.py output

args <- commandArgs(trailingOnly = TRUE)

if (args[1] == "-i") {
	# Load index data
	index_data = read.table("index_counts.txt", header=TRUE)

	# Correct python-to-R header parsing
	index_depth = ncol(index_data)/2
	colnames(index_data) <- c(-rev(seq(1,index_depth)), seq(1,index_depth))

	# Plot index data & save
	png("index_data.png")
	barplot(as.matrix(index_data), main="Base bias by position", 
		ylab= "UMI reads",
		xlab="Index", beside=TRUE, col=c("darkgreen","blue","red","yellow"))
	par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
	legend("topright", inset=c(-0.1,0), c(row.names(index_data)), cex=0.6, 
       bty="n", fill=c("darkgreen","blue","red","yellow"))
	dev.off()
}

# Read in pair analysis
pair_data = read.table("pair_counts.txt", header=TRUE)

# Correct python-to-R header parsing
pair_depth = ncol(pair_data)
colnames(pair_data) <- seq(1,pair_depth)

# Plot
png("pair_counts.png")
barplot(as.matrix(pair_data), main="Base bias by position", ylab= "UMI reads",
        xlab="Index", beside=TRUE, col=rainbow(16))
legend("topright", inset=c(-0.1,0), c(row.names(pair_data)), cex=0.6, 
       bty="n", fill=rainbow(16))
dev.off()

# Read in compressed pair analysis
c_pair_data = read.table("pair_counts_comp.txt", header=TRUE)

# Correct python-to-R header parsing
c_pair_depth = ncol(c_pair_data)
colnames(c_pair_data) <- seq(1,c_pair_depth)

# Plot
png("pair_counts_combined.png")
barplot(as.matrix(c_pair_data), main="Base bias by position; Combined", ylab= "UMI reads",
        xlab="Index", beside=TRUE, col=rainbow(10))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", inset=c(-0.1,0), c(row.names(c_pair_data)), cex=0.6, 
       bty="n", fill=rainbow(10))
dev.off()

# Compliment or antagonist base pairing
basepairs <- c("AT","TA","GC","CG")
non_basepairs <- c("AA","AG","AC", "TT", "TG", "TC", "GG","GA","GT" ,"CC","CT","CA")
identical <- c("AA","GG","TT","CC")
pair_data_compl <- pair_data[basepairs,]
pair_data_antg <- pair_data[non_basepairs,]
pair_data_same <- pair_data[identical,]

png("pair_compliments.png")
barplot(as.matrix(pair_data_compl),
        main="Paired Base bias by position\n(Complimentary)", 
        ylab= "UMI reads",
        xlab="Index", beside=FALSE, col=rainbow(4))
legend("topright", inset=c(-0.1,0), c(row.names(pair_data_compl)), cex=0.6, 
       bty="n", fill=rainbow(4))
dev.off()

png("pair_antags.png")
barplot(as.matrix(pair_data_antg),
        main="Paired Base bias by position\n(Antagonist)", 
        ylab= "UMI reads",
        xlab="Index", beside=FALSE, col=rainbow(12))
legend("topright", inset=c(-0.1,0), c(row.names(pair_data_antg)), cex=0.6, 
       bty="n", fill=rainbow(12))
dev.off()

png("pair_same.png")
barplot(as.matrix(pair_data_same),
        main="Paired Base bias by position\n(Identical)", 
        ylab= "UMI reads",
        xlab="Index", beside=FALSE, col=rainbow(4))
legend("topright", inset=c(-0.1,0), c(row.names(pair_data_same)), cex=0.6, 
       bty="n", fill=rainbow(4))
dev.off()

