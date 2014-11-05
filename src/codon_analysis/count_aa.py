import sys
from collections import defaultdict

filename = sys.argv[1]

aa_counts = defaultdict(int)
aa_site_counts = defaultdict(int)

for line in open(filename):
    line = line.strip().split("\t")

    for i in range(1,len(line)):
        aa = line[i]
        aa_counts[aa] += 1 
        aa_site_counts[str(i) + " " + aa] += 1

for k,v in aa_site_counts.iteritems():
    print k, v
