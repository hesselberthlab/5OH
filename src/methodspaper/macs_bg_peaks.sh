# yeast genome size
genomesize=12e6
CHROM_SIZES=~/ref/genomes/sacCer1/sacCer1.chrom.sizes

extsize=25
sample=SP8.assembly-sacCer1.align-uniq.strand.all.CPMs
bedgraph=${sample}.bg
exp_name=peaks/${sample}

# broadPeak autosql
ucscdir=/vol1/software/modules-sw/ucsc/build/v286
asfile=$ucscdir/kent/src/hg/lib/encode/broadPeak.as

min_peaklength=25

narrowpeak="${exp_name}.narrowPeak"
bigbed=${exp_name}_narrowpeak.bb

macs2 bdgpeakcall -i $bedgraph \
            -o $narrowpeak \
            --no-trackline \
            -l $min_peaklength

# sometimes the score exceeds the maximum (1000) defined by the
# broadPeak spec. Reformat the broadPeak file to covert > 1000
# to 1000
narrowpeak_tmpfile="$narrowpeak.tmp"
awk 'BEGIN {OFS="\t"} \
    { if ($5 > 1000) $5 = 1000; print $0}' \
    < $narrowpeak \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
    > $narrowpeak_tmpfile
mv $narrowpeak_tmpfile $narrowpeak

bedToBigBed -type=bed6+ -as=$asfile \
    $narrowpeak $CHROM_SIZES $bigbed

cp $bigbed ~/public_html/projects/5OH/methodspaper/${bigbed}
