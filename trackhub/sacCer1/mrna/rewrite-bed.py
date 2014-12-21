#! /usr/bin/env python

import sys
import ipdb

# load map
names = {}
for line in file(sys.argv[1]):
   
    sysname, common = [f.strip() for f in line.strip().split('\t')]
    names[sysname] = common

for line in file(sys.argv[2]):
    
    fields = line.strip().split('\t')

    name = fields[3]

    if name in names:
        common = names[name]

        fields[3] = common

    # write extra column for index
    fields.append(name)

    print '\t'.join(fields)
