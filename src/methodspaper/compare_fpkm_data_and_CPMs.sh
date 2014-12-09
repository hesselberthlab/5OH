<<DOC
Mapping pos/neg bedgraph data to FPKM file, ftw.
DOC

resultsdir=~/projects/5OH/results/methodspaper/
samplename=SP8.assembly-sacCer1.align-all
sampleprefix=${resultsdir}${samplename}
posbg=${sampleprefix}.strand.pos.CPMs.bg
negbg=${sampleprefix}.strand.neg.CPMs.bg
stat=sum
output=${sampleprefix}.FPKM_${stat}.tab

fpkmfile=/vol2/home/speach/ref/genomes/sacCer1/exp.fpkm.bed

#pos strand
awk '$6 == "+"' $fpkmfile \
    | bedtools map -a - -b $posbg -c 4 -o $stat \
    | awk '$7 != "."' \
    > $output

#neg strand
awk '$6 == "-"' $fpkmfile \
    | bedtools map -a - -b $negbg -c 4 -o $stat \
    | awk '$7 != "."' \
    >> $output

