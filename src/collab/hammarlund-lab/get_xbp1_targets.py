#! /usr/bin/env python

wormref = "wormRef.tab"
new_xbp1 = "new_xbp1_targets.tab"

for i,line in enumerate(open(new_xbp1)):
    if i == 1: continue
    gene_id = line.strip()
    printed = False
    for line in open(wormref):
        try:
            nm, gene, gene_id2 = line.strip().split("\t")
        except: continue
        if gene_id == gene_id2:
            print line.strip()
            printed = True
    if not printed:
        print "\t".join(["NA","NA",gene_id])
  