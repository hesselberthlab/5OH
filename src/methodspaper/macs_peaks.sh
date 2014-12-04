<<DOC
MACS peak-calling for 5OH-Seq data
DOC

# yeast genome size
genomesize=12e6
CHROM_SIZES=~/ref/genomes/sacCer1/sacCer1.chrom.sizes

sample=SP8.assembly-sacCer1.align-uniq
bam=${sample}.bam
exp_name=peaks/${sample}

# broadPeak autosql
ucscdir=/vol1/software/modules-sw/ucsc/build/v286
asfile=$ucscdir/kent/src/hg/lib/encode/broadPeak.as

macs2 callpeak -t $bam \
    -n $exp_name \
    --keep-dup all \
    --gsize $genomesize \
    --nomodel \
    --extsize 5 \
    -s 25 \
    --call-summits

peaktype=narrowPeak
peakfile=${exp_name}_peaks.${peaktype}
rawpeakfile=${exp_name}_peaks.raw.${peaktype}
bigbed=${exp_name}.ext5.s25.callsumits_${peaktype}.bb

cp $peakfile $rawpeakfile 

# sometimes the score exceeds the maximum (1000) defined by the
# broadPeak spec. Reformat the broadPeak file to covert > 1000
# to 1000
peak_tmpfile="$peakfile.tmp"
awk 'BEGIN {OFS="\t"} \
    { if ($5 > 1000) $5 = 1000; print $0}' \
    < $peakfile \
    | awk -v extsize=$extsize '$3 - $2 > extsize' \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
    > $peak_tmpfile
mv $peak_tmpfile $peakfile

bedToBigBed -type=bed6+4 -as=$asfile \
    $peakfile $CHROM_SIZES $bigbed

cp $bigbed ~/public_html/projects/5OH/methodspaper/${bigbed}
