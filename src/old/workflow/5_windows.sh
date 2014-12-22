# Creates an intersection table of normalized reads with mrna windows file
# annnnnd UTR windows file
MRNAWINDOWS="/vol2/home/speach/ref/genomes/sacCer1/sacCer1.mrna.20windows.bed"
UTRWINDOWS="/vol2/home/speach/ref/genomes/sacCer1/sacCer1.UTRs.2window.bed"
NORM_POSBGFILE="ATTGGC_S4.norm.pos.bg"
NORM_NEGBGFILE="ATTGGC_S4.norm.neg.bg"
WINDOW="ATTGGC_S4.window.full.tab"

# Section creating mRNA windows
bedtools intersect -a $NORM_POSBGFILE -b $MRNAWINDOWS -wao \
        | awk 'BEGIN{OFS="\t"; print "chr\tstart\tstop\tcount\tgene\tbin\tstrand\tcat"} $10=="+" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"+","exon"}' \
        > $WINDOW

bedtools intersect -a $NORM_NEGBGFILE -b $MRNAWINDOWS -wao \
        | awk 'BEGIN{OFS="\t"} $10=="-" {split($8, a, "_"); print $1, $2, $3, $4, a[1], a[2],"-","exon"}' \
        >> $WINDOW

# Section creating UTR windows
bedtools intersect -a $NORM_POSBGFILE -b $UTRWINDOWS -wao \
    | awk 'BEGIN{OFS="\t"} $10=="+" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"+",a[2]}' \
    >> $WINDOW

bedtools intersect -a $NORM_NEGBGFILE -b $UTRWINDOWS -wao \
    | awk 'BEGIN{OFS="\t"} $10=="-" {split($8,a,"_"); print $1, $2, $3, $4, a[1],a[3],"-",a[2]}' \
    >> $WINDOW
    
