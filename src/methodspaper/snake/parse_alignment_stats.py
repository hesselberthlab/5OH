'''
Parse bowtie alignment output data generated via snakemake
'''

import os
import sys

results_dir = sys.argv[1]

# Print header for information
print "\t".join(["sample", "assembly", "alignment", "processed",
"aligned", "failed",
"suppressed"])


# Cycle through files and print out 
for filename in os.listdir(results_dir):
    filefields = filename.strip().split(".")
    if len(filefields) == 5 and filefields[3] == "alignstats":
        sample, assembly, alignment = filefields[:3]
        assembly = assembly.split("-")[1]
        alignment = alignment.split("-")[1]

        read_data = [sample, assembly, alignment]

        for line in open(filename):
            if line.startswith("Reported "): break # final line of output
            read_datum = line.strip().split(": ")[1].split(" (")[0]
            read_data.append(read_datum)

        if len(read_data) != 7:
            read_data.append("NA")
        print "\t".join(read_data)



