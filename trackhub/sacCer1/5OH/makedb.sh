#! /usr/bin/env bash

set -o nounset -o pipefail -o errexit -x

mainstanza="
track 5OH-seq\n
compositeTrack on\n
shortLabel 5OH-seq\n
longLabel 5OH-seq Coverage and Peaks tracks\n
subGroup1 view Views COV=Coverage\n
subGroup2 strain Strain WT=wild-type delXrn1=delXrn1\n
subGroup3 treatment Treatment NONE=None Tm=Tunicamycin SAP=phosphatase DMSO=DMSO DMSO.SAP=DMSO+phosphatase Tm.SAP=Tunicamycin+phosphatase PreFrag=Pre-fragmentation\n
subGroup4 strand Strand pos=pos neg=neg all=all\n
subGroup5 replicate Replicate NONE=none Rep1=Rep1 Rep2=Rep2\n
subGroup6 signal Signal CPMs=CPMs counts=raw-counts\n
dimensions dimX=treatment dimY=strain dimA=strand dimB=replicate dimC=signal\n
filterComposite dimA dimB dimC\n
sortOrder view=+ strand=+\n
type bed 6 .\n
\n
\ttrack 5OH-Coverage\n
\tview COV\n
\tparent 5OH-seq\n
\tshortLabel 5OH-seq Coverage\n
\tvisibility full\n
\ttype bigWig\n
\tautoScale on\n"

echo -e $mainstanza

# XXX Yo SP: What is "align.std"???

strains=(WT delXrn1)
treatments=(None PreFrag Tm Tm.SAP DMSO DMSO.SAP)
replicates=(Rep1 Rep2)
strands=(pos neg all)
signals=(counts CPMs)

for strand in ${strands[@]}; do
  for strain in ${strains[@]}; do
    for treatment in ${treatments[@]}; do
      for signal in ${signals[@]}; do
        for rep in ${replicates[@]}; do

        fname="$strain.$treatment.$rep.align.std.strand.$strand.$signals.bw"

        if [[ -f $fname ]]; then

            if [[ $strand == "pos" ]]; then
                color="202,0,32"
            else
                color="5,113,176"
            fi

            stanza="
                    \t\ttrack $strain.$treatment.$rep.$strand.$signal\n
                    \t\tparent 5OH-Coverage\n
                    \t\tsubGroups view=COV strain=$strain strand=$strand treatment=$treatment replicate=$rep signal=$signal\n
                    \t\tbigDataUrl 5OH/$fname\n
                    \t\tshortLabel $strain.$treatment.$rep\n
                    \t\tlongLabel 5OH-seq coverage strain=$strain treatment=$treatment strand=$strand replicate=$rep\n
                    \t\tmaxHeightPixels 30:30:10\n
                    \t\tcolor $color\n
                    \t\ttype bigWig\n
            "
            echo -e $stanza
        # else
            # echo ">> file $fname not found"
        fi

        done
      done            
    done
  done
done
