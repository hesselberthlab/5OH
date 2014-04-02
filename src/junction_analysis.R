# Creates simple R plots from junction_analysis.py output

# Load index data
index_data = read.table("index_counts.txt", header=TRUE)

# Correct python-to-R header parsing
index_depth = ncol(index_data)/2
colnames(index_data) <- c(-rev(seq(1,index_depth)), seq(1,index_depth))

# Plot index data
barplot(as.matrix(index_data), main="Base bias by position", ylab= "UMI reads",
        xlab="Index", beside=TRUE, col=c("darkgreen","blue","red","yellow"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", inset=c(-0.1,0), c(row.names(index_data)), cex=0.6, 
       bty="n", fill=c("darkgreen","blue","red","yellow"))

# Read in pair analysis
pair_data = read.table("pair_counts.txt", header=TRUE)

# Correct python-to-R header parsing
pair_depth = ncol(pair_data)
colnames(pair_data) <- seq(1,pair_depth)

# Plot
barplot(as.matrix(pair_data), main="Base bias by position", ylab= "UMI reads",
        xlab="Index", beside=TRUE, col=rainbow(16))
legend("topright", inset=c(-0.1,0), c(row.names(pair_data)), cex=0.6, 
       bty="n", fill=rainbow(16))

# Read in compressed pair analysis
c_pair_data = read.table("pair_counts_comp.txt", header=TRUE)

# Correct python-to-R header parsing
c_pair_depth = ncol(c_pair_data)
colnames(c_pair_data) <- seq(1,c_pair_depth)

# Plot
barplot(as.matrix(c_pair_data), main="Base bias by position; Combined", ylab= "UMI reads",
        xlab="Index", beside=TRUE, col=rainbow(10))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", inset=c(-0.1,0), c(row.names(c_pair_data)), cex=0.6, 
       bty="n", fill=rainbow(10))

