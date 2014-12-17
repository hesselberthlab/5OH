import sys
from subprocess import Popen, PIPE


colnum = int(sys.argv[1]) -1
inputfile = sys.argv[2]
namemap = "/vol2/home/speach/ref/genomes/sacCer1/sacCer1.sysname-genename-map.txt"

for line in open(inputfile):
    line = line.strip().split("\t")
    name = line[colnum]

    bashcommand = "%s %s %s" % ("grep", name, namemap)
    process = Popen(bashcommand.split(), stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    
    if len(stdout) == 0:
        newname = name
    else:
        newname = stdout.strip().split("\t")[-1]
   
    line[colnum] = newname

    print "\t".join(line)
