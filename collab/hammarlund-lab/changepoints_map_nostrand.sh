#! /usr/bin/env bash

<<DOC
1) Generate bedgraphs with ALL bases in genome
2) Map onto reference genome
3) Caluclate changepoints

input: DATA + "{sample}.tophat_/accepted_hits.bam"
params: REFBED
output: RESULT + "{sample}.exonmappings.tab"
DOC

# Read Arguments
bamfile=$1
refbed=$2
mapfile=$3

sample=$(basename $mapfile .exonmappings.tab)
results=$(dirname $mapfile)
allbasetable="${results}/${sample}.allbase.tab"
allbasebg="${results}/${sample}.allbase.bg"

# Report coverage at every base in the genome, minus introns
# Output is 3-columns: chrom, base, coverage
bedtools genomecov -ibam $bamfile -d > $allbasetable

# Convert to bedgraph format by adding in "stop" column
awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $allbasetable > $allbasebg

# Map to genes
awk 'BEGIN{OFS="\t"} $6 == "+"' $refbed \
    | cut -f 1-6 \
    | bedtools map -a - -b $allbasebg -c 4 -o collapse \
    > $mapfile
    
awk 'BEGIN{OFS="\t"} $6 == "-"' $refbed \
    | cut -f 1-6 \
    | bedtools map -a - -b $allbasebg -c 4 -o collapse \
    >> $mapfile