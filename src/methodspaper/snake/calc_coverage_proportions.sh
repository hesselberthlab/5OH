# Scripts for converting an "all-base bedgraph" where the score is the number of CPMs (or raw counts)
# to an "all-base bedgraph" where the score is the proportion of total reads within a gene model
# that occur at a given base...
# Simple implementation doesn't account for overlapping gene models (TODO)

devel=~/devel/5OH/src/codon_analysis/
CHROMSIZES=~/ref/genomes/sacCer1/sacCer1.chrom.sizes

# generate allbase bg (TODO: compare to previous mappings)
sample=$1
subdir=proportions/

bam=$1
genebed=$2
subdir=$3
propsbed=$4

sample=$(basename $bam .bam)

# Create all these lovely stranded files of beauty
allbasetable=${subdir}${sample}.allbase.tab
allpostable=${subdir}${sample}.allpos.tab
allnegtable=${subdir}${sample}.allneg.tab

allbasebg=${subdir}${sample}.allbase.bg
allposbg=${subdir}${sample}.allpos.bg
allnegbg=${subdir}${sample}.allneg.bg

if [[ ! -d $subdir ]] ; then
    mkdir $subdir
fi

bedtools genomecov -5 -ibam $bam -d > $allbasetable
bedtools genomecov -5 -ibam $bam -d -strand + > $allpostable
bedtools genomecov -5 -ibam $bam -d -strand - > $allnegtable

awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $allbasetable > $allbasebg
awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $allpostable > $allposbg
awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $3}' $allnegtable > $allnegbg

# map bg to genebed to get coverage on gene model
# do it stranded, because life is grand
#posgenemapfile=${subdir}${sample}.genemap.pos.tab
#neggenemapfile=${subdir}${sample}.genemap.neg.tab
genemapfile=${subdir}${sample}.genemap.tab

awk '$6 == "+"' $genebed \
    | bedtools map -a - -b $allposbg -c 4 -o collapse > $genemapfile

awk '$6 == "-"' $genebed \
    | bedtools map -a - -b $allnegbg -c 4 -o collapse >> $genemapfile

python ${devel}calc_coverage_proportions.py $genemapfile > ${propsbed}.unsorted
sort -k1,1 -k2,2n ${propsbed}.unsorted > $propsbed
