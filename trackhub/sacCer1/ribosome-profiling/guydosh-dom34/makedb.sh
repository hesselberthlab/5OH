#! /usr/bin/env bash

mainstanza="
track Ribosome-Footprints\n
compositeTrack on\n
shortLabel Ribosome footprints\n
longLabel Ribosome Footprints Coverage and Peaks tracks\n
subGroup1 view Views COV=Coverage\n
subGroup2 strain Strain WT=wild-type DOM=dom34KO SKI=ski2KO HBS=hbs1KO\n
subGroup3 treatment Treatment NONE=None CHX=CHX MG=high-Mg2+ DIAM=diamde GP=GMP-PNP CGP=CHX+GMP-PNP LG=low-glucose TR=suppressor-tRNA AT=3-AT\n
subGroup4 strand Strand pos=pos neg=neg\n
subGroup5 size Size UNF=unfractionated SHR=short DIS=disome\n
dimensions dimX=treatment dimY=strain dimA=strand dimB=size\n
filterComposite dimA dimB\n
sortOrder view=-\n
type bed 6 .\n
\n
\ttrack Coverage\n
\tview COV\n
\tparent Ribosome-Footprints\n
\tshortLabel Coverage\n
\tvisibility full\n
\ttype bigWig\n
\tautoScale on\n"

echo -e $mainstanza

accessions=$(cut -f1 accessions.txt)
strands=(pos neg)

for acc in $accessions; do

    descrip=$(grep $acc accessions.txt | cut -f2)
    strain=$(grep $acc accessions.txt | cut -f3)
    treatment=$(grep $acc accessions.txt | cut -f4)
    size=$(grep $acc accessions.txt | cut -f5)

    for strand in ${strands[@]}; do

        if [[ $strand == "pos" ]]; then
            color="202,0,32"
        else
            color="5,113,176"
        fi

        stanza="
                \t\ttrack $acc.$strand.coverage\n
                \t\tparent Coverage\n
                \t\tsubGroups view=COV strain=$strain strand=$strand treatment=$treatment size=$size\n
                \t\tbigDataUrl ribosome-profiling/guydosh-dom34/$acc.$strand.bw\n
                \t\tshortLabel $strain.$treatment.$size.$strand\n
                \t\tlongLabel $acc $descrip $strand\n
                \t\tmaxHeightPixels 30:30:10\n
                \t\tcolor $color\n
                \t\ttype bigWig\n
        "
        echo -e $stanza
    done
done
