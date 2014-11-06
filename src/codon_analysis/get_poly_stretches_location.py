#! /usr/bin/env python

'''get_poly_stretches_location.py: find occurences of poly-aa stretches in gene model
    output: gene and distance to start of poly-aa stretch;
            multiple poly-aa stretches are comma-separated in 2nd column'''
            
# this code is pretty ugly, but whatever.

from subprocess import Popen, PIPE
#import ipdb


aminos = ["E"]
#aminos = ["D", "E", "[DE]", "R", "K", "[RK]"]
protein_list = '/vol2/home/speach/ref/genomes/sacCer1/orf_trans.tab'
full_genome_file = '/vol2/home/speach/ref/genomes/sacCer1/sacCer1.sgdGene.tab'
POLY_LENGTH = 5 # length of amino acid expansions
SPAN = 50 # length to search upstream and downstream of poly-aa stretch

def get_location(aminos, protein_list, mismatch=0):
    # set shell command call; agrep allows for a regex with a predifined number of mismatches
    if mismatch == 0:
        grep = "grep"
    else:
        grep = "agrep -%s" % mismatch
    
    # loop through amino acids, grep genes via bash subprocess
    for i in range (1,5):
        for amino in aminos:
            polyaa = amino * POLY_LENGTH * i
            polyaa_name = amino + str(POLY_LENGTH * i)
            
            bashcommand = "%s %s %s | wc" % (grep, polyaa, protein_list)
            process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
            genes = stdout.split("\n")[:-1]
            
            for gene in genes:
                gene, seq = gene.split("tab:")[-1].split("\t\t")
                poly_start = []
                
                # unique poly-aa stretch
                if seq.count(polyaa) == 1:
                    poly_start.append(str(seq.index(polyaa)*3))
                    
                # if there is more than one occurance of poly-aa stretch,
                # append to list if it is a unique region.
                # ie, two 5aa E stretches that are part of a continuous
                # 10aa E stretch would only be counted once
                else:
                    prev_end = 0
                    for i in range(0,seq.count(polyaa)):
                        start = seq.index(polyaa, prev_end)
                        if start == prev_end:
                            prev_end = start + POLY_LENGTH
                            continue
                        poly_start.append(str(start*3))
                        prev_end = start + POLY_LENGTH
                
                for start in poly_start:
                    get_genome_coords(gene,polyaa_name,start,SPAN)
                #print "\t".join([gene, ",".join(poly_start)])

def get_genome_coords(gene,polyaa_name,polyaa_start,span):
    bashcommand = "grep %s %s" % (gene, full_genome_file)
    process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    gene, chrom, strand, gene_start, gene_stop = stdout.strip().split("\t")[:5]
    gene_start = int(gene_start)
    gene_stop = int(gene_stop)
    polyaa_start = int(polyaa_start)
    
    if strand == "+":
        query_stop = gene_start + polyaa_start + span + 1
        query_start = gene_start + polyaa_start - span
    
    if strand == "-":
        query_start = gene_stop - polyaa_start - span - 1
        query_stop = gene_stop - polyaa_start + span
    
    bedline = [chrom, str(query_start), str(query_stop), gene, polyaa_name, strand]
    print "\t".join(bedline)


get_location(aminos, protein_list, 0)
    
