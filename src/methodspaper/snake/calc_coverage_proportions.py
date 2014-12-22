'''Calculate coverage proportions from 5OH-Seq data mapped to all bases of mRNA gene models'''

import sys

genemapfile = sys.argv[1]

for line in open(genemapfile):
    # parse line
    chrom, gene_start, gene_stop, gene, score, strand, counts = line.strip().split("\t")
    gene_start = int(gene_start)
    gene_stop = int(gene_stop)
    
    # Split map counts
    countarray = counts.split(",")
    if countarray == ["."]: continue # no data maps here
    else:
        countarray = map(int, countarray)
        if max(countarray) < 15: continue # low counts, eff that
    
    total_counts = float(sum(countarray))
    
    for i in range(gene_stop-gene_start):
        count = countarray[i]
        
        if count == 0: # ain't nothing here
            continue
        
        else:
            start = gene_start + i
            stop = start + 1
            proportion = count / total_counts
            print "\t".join([chrom, str(start), str(stop), gene, str(proportion), strand])
    
