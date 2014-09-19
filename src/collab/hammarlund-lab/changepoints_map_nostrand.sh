#! /usr/bin/env bash

<<DOC
1) Generate bedgraphs with ALL bases in genome
2) Map onto reference genome
3) Caluclate changepoints

input: DATA + "{sample}.tophat_/accepted_hits.bam"
params: REFBED
output: RESULT + "{sample}.pos.mappings.tab", RESULT + "{sample}.neg.mappings.tab" 
DOC

# Read Arguments
bamfile=$1
refbed=$2
posmapfile=$3
negmapfile=$4

sample=$(basename $posmapfile .pos.mappings.tab)
results=$(dirname $posmapfile)
allbasetable="${results}/${sample}.allbase.tab"
allbasebg="${results}/${sample}.allbase.bg"

# Report coverage at every base in the genome
# Output is 3-columns: chrom, base, coverage
bedtools genomecov -ibam $bamfile -d > $allbasetable

# TODO: Fix strand issue?
#bedtools genomecov -ibam $bamfile -d -strand + > $posallbasetable
#bedtools genomecov -ibam $bamfile -d -strand - > $negallbasetable

# Convert to bedgraph format by adding in "stop" column
awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $allbasetable > $allbasebg

# TODO: check +/- 1 errors with strandedness.
#awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $posallbasetable > $posallbasebg
#awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $negallbasetable > $negallbasebg

# Map to genes
awk 'BEGIN{OFS="\t"} $6 == "+"' $refbed \
    | cut -f 1-6 \
    | bedtools map -a - -b $allbasebg -c 4 -o collapse \
    > $posmapfile
    
awk 'BEGIN{OFS="\t"} $6 == "-"' $refbed \
    | cut -f 1-6 \
    | bedtools map -a - -b $allbasebg -c 4 -o collapse \
    > $negmapfile