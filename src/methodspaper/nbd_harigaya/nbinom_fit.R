#! /usr/bin/env Rscript

# nbinom_fit.R
# usage: Rscript nbinom_fit.R [infile] [outfile]

library(MASS)
args=commandArgs(TRUE)
infile=as.character(args[1])
outfile=as.character(args[2])

input=readLines(infile)
n=length(input)
# Create a matrix where the p_values will be stored.
prob=matrix(0, 0, 2)
# Create a vector where names of fitted genes will be stored.
fit=c()
# Create a vector where names of unfitted genes will be stored.
unfit=c()

for (i in 1:n)
{
# Split each line into a list and convert the list into a character vector.
	line=unlist(strsplit(input[i], "\t"))
# Store the gene name and profile separately.
	name=line[1]
	profile=as.integer(line[-1])
# Fit the profile to a negative binomial distribution.
	model=try(fitdistr(profile, "Negative Binomial"), silent = TRUE)
	if(attr(model,"class") == "fitdistr")
	{
# Add the gene name to the fit vector.
		fit=c(fit, name)
# Extract the model information, size and mu.
		parameters=coef(model)
# Compute p_values at each site and store in a numeric vector.
		p_values=dnbinom(profile, size=parameters[1], mu=parameters[2])
# Make a matrix.
		m=length(p_values)
		p_mat=matrix(0, m, 2)
		for(j in 1:m)
		{
			p_mat[j,1]=paste(name,"_",j,sep="")
			p_mat[j,2]=p_values[j]
		}
# Add the matrix to the prob matrix.
		prob=rbind(prob, p_mat)
	}else
	{
# Add the gene name to the unfit vector.
		unfit=c(unfit, name)		
	}
}

output=file(outfile, "w")
write.table(prob, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", eol="\n")
close(output)

outfile2=paste(outfile,".fit",sep="")
output=file(outfile2, "w")
write.table(fit, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", eol="\n")
close(output)

outfile3=paste(outfile,".unfit", sep="")
output=file(outfile3, "w")
write.table(unfit, file=output, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", eol="\n")
close(output)

q()
	
	





