<<DOC
Intersect MACS narrow peak file with proportions bed to get an
R-compatible table of:
MACS peak: chrom, start, stop, gene reads
max proportion site: proportion, chrom, position, strand
DOC

sample=SP8.assembly-sacCer1.align-uniq
proportionsbed=${sample}.maxprops.bed
narrowpeak=${sample}_peaks.raw.narrowPeak
output=${sample}_peaks.intersect.tab

bedtools intersect -a $proportionsbed -b $narrowpeak
    | awk 'BEGIN{OFS="\t"} $7 != "." {print $7, $8, $9, $4, $11, $5, $1, $2, $6}' \
    > $output
