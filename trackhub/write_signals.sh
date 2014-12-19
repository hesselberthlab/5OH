#! /usr/bin/env bash

# write signal track section for hub
# see this file for reference:
# /vol4/home/brownj/projects/polya/results/common/hub/hg19/trackDb.txt

mainstanza="
track 5OH-seq
compositeTrack on
shortLabel Coverage & Peaks
longLable Coverage and Peaks tracks
subGroup1 view Views COV=Coverage PKS=Peaks
subGroupXXX strain Strains WT=WT XRN1=delXrn1
subGroup2 treatment Treatment XXX
subGroupXXX rep Replicate REP1=Rep1 REP2=Rep2
subGroupXXX strand Strand POS=pos NEG=neg
dimensions dimX=strain dimY=treatment dimA=rep dimB=strand
filterComposite dimA dimB
sortOrder view=-
type bed 6 .
"

strains=(WT delXrn1)
treatments=(DMSO DMSO.SAP Tm Tm.SAP)
signals=(counts CPMs)
strands=(pos neg)
reps=(Rep1 Rep2)

# write out stanzas
for strain in ${strains[@]}; do
  for treatment in ${treatments[@]}; do
    for signal in ${signals[@]}; do
      for strand in ${strands[@]}; do
        for rep in ${reps[@]}; do
            
        done
      done
    done
  done
done
