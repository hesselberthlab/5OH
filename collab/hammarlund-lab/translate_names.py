#! /usr/bin/env python

# Prints systematic name of c. elegans genes

import re

targets_filename = "known_xbp1_targets.tab"
dictionary_filename = "refseq.wormgene.name.map.tab"
genediff_filename = "gene_exp.diff"

xbp1_targets = []

for line in open(targets_filename):
    # most genes not in dictionary as proper isoform
    # can look at all versions by stripping dot notation
    gene = line.strip() #.split(".")[0]
    
    for line in open(dictionary_filename):
        line = line.strip()
        if line.endswith(gene):
            #print line
            systematic_name = line.split("\t")[0]
            xbp1_targets.append(systematic_name)

for gene in xbp1_targets:
    for line in open(genediff_filename):
        if re.search(gene,line):
            print line.strip()


