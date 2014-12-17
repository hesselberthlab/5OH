#! /usr/bin/env python

'''get_poly_stretches_location.py: find occurences of poly-aa stretches
    output: bedfile of some number of bp upstream and downstream of all poly-aa stretches'''

from subprocess import Popen, PIPE
import sys
import re # reg is my fav
import regex #regex is my FAV fav

amino = sys.argv[1]
aminos = []
aminos.append(amino)

if amino == "-all":
    aminos_string = "ACDEFGHIKLMNPQRSTVWY"
    aminos = list(aminos_string)
    aminos.append("[DE]")
    aminos.append("[RK]")
    aminos.append("[RKH]")
    aminos.append("[ST]")
    aminos.append("[NQ]")
    
#aminos = ["D", "E", "[DE]", "R", "K", "[RK]"]
protein_list = '/vol2/home/speach/ref/genomes/sacCer1/orf_trans.tab'
full_genome_file = '/vol2/home/speach/ref/genomes/sacCer1/sacCer1.sgdGene.tab'
POLY_LENGTH = 1 # length of amino acid expansions
SPAN = 75 # length to search upstream and downstream of poly-aa stretch
mismatch = 5

def get_location(aminos, protein_list, mismatch=0):
    # set shell command call; agrep allows for a regex with a predifined number of mismatches
    if mismatch == 0:
        grep = "grep"
    else:
        grep = "agrep -%s" % mismatch
    
    # loop through amino acids, grep genes via bash subprocess
    for i in range (10,11):
        for amino in aminos:
            polyaa = amino * POLY_LENGTH * i
            polyaa_name = amino + ";ML:" + str(POLY_LENGTH * i) + ";MM:" + str(mismatch)
    
    #for i in range (1,5):
    #    for amino in aminos:
    #        polyaa = amino * POLY_LENGTH * i
    #        polyaa_name = amino + str(POLY_LENGTH * i)
            
            bashcommand = "%s %s %s" % (grep, polyaa, protein_list)
            process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
            genes = stdout.split("\n")[:-1]
            
            for gene in genes:
                gene, seq = gene.split("\t\t")
                
                poly_start = identify_polystart(polyaa, seq, mismatch)
                
                for start in poly_start:
                    get_genome_coords(gene,polyaa_name,start,SPAN)
                #print "\t".join([gene, ",".join(poly_start)])

def identify_polystart(polyaa, seq, mismatch): 
    '''Identify start of poly-aa stretch'''
    
    # create regex search pattern
    pattern = '(%s){s<=%s}' % (polyaa, mismatch)

    # count occurances in sequence.
    # bestmatch flag causes regex to find matches closest to the pattern,
    # not just the first occurence
    patterns_found = regex.findall(pattern, seq, regex.BESTMATCH)


    poly_start = []
    
    # if there is a unique polyaa stretch in sequence
    if len(patterns_found) == 1:
          motif_pattern = patterns_found.pop(0)
          poly_start.append(seq.index(motif_pattern)*3)
                    
    # if there is more than one occurance of poly-aa stretch,
    # append to list if it is a unique region.
    # ie, two 5aa E stretches that are part of a continuous
    # 10aa E stretch would only be counted once
    else:
          prev_end = 0
          for j in range(0,len(patterns_found)):
                # remove patterns from front of list
                motif_pattern = patterns_found.pop(0)

                start = seq.index(motif_pattern, prev_end)
                if start == prev_end:
                    prev_end = start + POLY_LENGTH
                    continue
                poly_start.append(start*3)
                prev_end = start + POLY_LENGTH
    
    return poly_start


def get_genome_coords(gene,polyaa_name,polyaa_start,span):
    '''print bedline of span upstream and downstream of start of polyaa motif'''
    
    bashcommand = "grep %s %s" % (gene, full_genome_file)
    process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    if stdout:
        gene, chrom, strand, gene_start, gene_stop = stdout.strip().split("\t")[:5]
        gene_start = int(gene_start)
        gene_stop = int(gene_stop)
    else: return
    
    if strand == "+":
        query_stop = gene_start + polyaa_start + span + 1
        query_start = gene_start + polyaa_start - span
    
    if strand == "-":
        query_start = gene_stop - polyaa_start - span - 1
        query_stop = gene_stop - polyaa_start + span
    
    bedline = [chrom, str(query_start), str(query_stop), gene, polyaa_name, strand]
    print "\t".join(bedline)


# Main
for m in range(0, mismatch+1):
    get_location(aminos, protein_list, m)
