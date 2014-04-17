#! /usr/bin/env python

# Grab output statistics from alignment error files and umi-trimmed bedgraphs

import sys
from subprocess import Popen, PIPE

def get_alignment_stats(filename):
	for line in open(filename):
		if line.startswith('+ BAMPREFIX'):
			sample = line.strip().split("=")[1]
		if line.startswith('#'):
			line = line.strip()
			
			if 'processed' in line:
				processed = line.split(": ")[1]
			elif 'reported alignment' in line:
				aligned = line.split(": ")[1].split(" ")[0]
			elif 'failed' in line:
				failed = line.split(": ")[1].split(" ")[0]
			elif 'suppressed' in line:
				suppressed = line.split(": ")[1].split(" ")[0]
		
		umi_reads,sites = get_umi_reads_and_sites(sample)
	
	print "\t".join([sample, processed, aligned, failed, 
	    suppressed,umi_reads,sites])

def get_umi_reads_and_sites(sample):
    S_number = sample.split("_")[1]
    command = "awk '{sum+=$4} END{print sum}' *%s.bg" % (S_number)
    umi_reads = Popen(command, stdout=PIPE, shell=True).stdout.read()
    command = "wc -l *%s.bg" % (S_number)
    sites = Popen(command, stdout=PIPE, shell=True).stdout.read()
    
    return umi_reads.strip(), sites.strip().split(" ")[0]

	
def main():
	print "\t".join(["sample","processed","aligned","failed","suppressed","umireads","sites"])
	for filename in Popen("ls err-out/*.err", stdout=PIPE, shell=True).stdout:
		get_alignment_stats(filename.strip())

if __name__ == '__main__':
    sys.exit(main()) 
