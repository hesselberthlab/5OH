#! /usr/bin/env python

'''get_poly_stretches_location.py: find occurences of poly-aa stretches
    output: bedfile of some number of bp upstream and downstream of all poly-aa stretches'''

from subprocess import Popen, PIPE
import sys
import re #regex is my fav

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
    aminos.append("[RKED]")
    aminos.append("[RD]")
    aminos.append("[KD]")
    aminos.append("[KE]")
    aminos.append("[RE]")
    aminos.append("[RKE]")
    aminos.append("[RKD]")
    
#aminos = ["D", "E", "[DE]", "R", "K", "[RK]"]
protein_list = '/vol2/home/speach/ref/genomes/sacCer1/orf_trans.tab'
full_genome_file = '/vol2/home/speach/ref/genomes/sacCer1/sacCer1.sgdGene.tab'
POLY_INCREMENTS = 5 # length of amino acid expansions
SPAN = 75 # length to search upstream and downstream of poly-aa stretch

def get_location(aminos, protein_list, mismatch=0):
    # set shell command call; agrep allows for a regex with a predifined number of mismatches
    if mismatch == 0:
        grep = "grep"
    else:
        grep = "agrep -%s" % mismatch
    
    # loop through amino acids, grep genes via bash subprocess
    for motif_length in range (5,21,POLY_INCREMENTS):
        for amino in aminos:
            polyaa = amino * motif_length
            polyaa_name = amino + str(motif_length)
    
            bashcommand = "%s %s %s" % (grep, polyaa, protein_list)
            process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
            genes = stdout.split("\n")[:-1]
            
            for gene in genes:
                gene, seq = gene.split("\t\t")
                
                if polyaa.startswith("["):
                    poly_start = regex_identify_polystart(polyaa, seq,
                                                           motif_length)
                else:
                    poly_start = identify_polystart(polyaa, seq,
                                                    motif_length)
                
                for start in poly_start:
                    get_genome_coords(gene,polyaa_name,start,SPAN)
                #print "\t".join([gene, ",".join(poly_start)])

def identify_polystart(polyaa, seq, motif_length): 
    '''Identify start of poly-aa stretch'''
        
    poly_start = []
    
    # if there is a unique polyaa stretch in sequence
    if seq.count(polyaa) == 1:
          poly_start.append(seq.index(polyaa)*3)
                    
    # if there is more than one occurance of poly-aa stretch,
    # append to list if it is a unique region.
    # ie, two 5aa E stretches that are part of a continuous
    # 10aa E stretch would only be counted once
    else:
          prev_end = 0
          for j in range(0,seq.count(polyaa)):
                start = seq.index(polyaa, prev_end)
                if start == prev_end:
                    prev_end = start + motif_length
                    continue
                poly_start.append(start*3)
                prev_end = start + motif_length
    
    return poly_start

def regex_identify_polystart(polyaa, seq, motif_length): 
    '''Identify start of poly-aa stretch when regex poly-aa is used (ie, polyacidic, polybasic)'''
        
    poly_start = []
    
    allmatches = re.findall(polyaa,seq) # returns list of polyaa stretches
                                        # matching regex
    
    # if there is a unique polyaa stretch in sequence
    if len(allmatches) == 1:
          polyaa = allmatches[0]
          poly_start.append(seq.index(polyaa)*3)
                    
    # if there is more than one occurance of poly-aa stretch,
    # append to list if it is a unique region.
    # ie, two 5aa DE stretches that are part of a continuous
    # 10aa DE stretch would only be counted once
    else:
          prev_end = 0
          for i in range(0,len(allmatches)):
                start = seq.index(allmatches[i], prev_end)
                if start == prev_end:
                    prev_end = start + motif_length
                    continue
                poly_start.append(start*3)
                prev_end = start + motif_length
    
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

# make main
get_location(aminos, protein_list, 0)
    
