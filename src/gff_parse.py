#! /usr/bin/env python

# Parse gff file trimmed of header and FASTA

import sys
import re
from collections import defaultdict

gene_dict = defaultdict(str)

def gff_parse(filename):
    prev_start = 0
    prev2_start = 0     # parses out extra mRNA entry after intron
    prev_stop = 0
    
    for line in open(filename):
        fields = line.strip().split("\t")
        chrm, sgd, cat, start, stop, d1, strand, d2, descr = fields[:]
        descr = reg_parse(descr)
        chrm = edit_chrm(chrm)
        
        if cat == "region" or cat == "chromosome": continue     # nix crap
        #if start == prev_start or stop == prev_stop: continue  # prevent rep
        #if start == prev2_start: continue
        
        # grab understandable names
        if "gene=" in descr:
            ID = descr.split("ID=")[1].split(";")[0]
            name = descr.split("gene=")[1].split(";")[0]
            gene_dict[ID] = name
            
        elif "Name=" in descr:
            name = descr.split("Name=")[1].split(";")[0]
            name = name.split("_")[0]
            if name in gene_dict.keys():
                name = gene_dict[name]
        
        # add additional CDS classification
        if cat == "CDS" and "orf_classification=" in descr:
            cat+= "_" + descr.split("orf_classification=")[1]
        
        # change mito chrom to proper name
        if chrm == "chrmt":
            chrm = "chrM"
        
        # grab notes in long form
        #if "Note=" in descr:
        #    cat+= "\t" + descr.split(";")[3].split("=")[1]
            
        
        print "\t".join([chrm, start, stop, ":".join([cat,name]), "1", strand])
       
        prev2_start = prev_start
        prev_start = start
        prev_stop = stop
        
"""
def trim_duplicates(gff_parse_iterator):
    prev_coords = (0,0)
    prev2_coords = (0,0)
    
    for chrm,start,stop,cat,num,strand in gff_parse_iterator:
        if prev_cords == (start,stop) && prev2_coords == (start,stop):
            pass
"""   

def reg_parse(line):
    """ Parse "%n" characters back to text
        Params:
            line (str) - unedited
        Returns:
            line (str) - regex edited
    """
    line = re.sub(r'\%28','(',line)
    line = re.sub(r'\%29',')',line)
    line = re.sub(r'\%20',' ',line)
    
    return line
    
def edit_chrm(chrm):
    roman = chrm[3:]
    chrom_dict = {"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6,"VII":7,"VIII":8,
                    "IX":9,"X":10,"XI":11,"XII":12,"XIII":13,"XIV":14,
                    "XV":15,"XVI":16,"mt":"M"}
    digit = chrom_dict[roman]
    return "chr" + str(digit)
    

def main(args=sys.argv[1:]):
    filename = args[0]
    return gff_parse(filename)

if __name__ == '__main__':
    sys.exit(main())
