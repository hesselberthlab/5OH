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
posbasetable="${results}/${sample}.posbase.tab"
posbasebg="${results}/${sample}.posbase.bg"
negbasetable="${results}/${sample}.negbase.tab"
negbasebg="${results}/${sample}.negbase.bg"

# Report coverage at every base in the genome
# Output is 3-columns: chrom, base, coverage
bedtools genomecov -ibam $bamfile -d -strand + > $posbasetable
bedtools genomecov -ibam $bamfile -d -strand - > $negbasetable

# Convert to bedgraph format by adding in "stop" column
awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $posbasetable > $posbasebg
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' $negbasetable > $negbasebg

# Map to genes
awk 'BEGIN{OFS="\t"} $6 == "+"' $refbed \
    | cut -f 1-6 \
    | bedtools map -a - -b $posbasebg -c 4 -o collapse \
    > $mapfile
    
awk 'BEGIN{OFS="\t"} $6 == "-"' $refbed \
    | cut -f 1-6 \
    | bedtools map -a - -b $negbasebg -c 4 -o collapse \
    >> $mapfile