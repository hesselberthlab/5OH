#! /usr/bin/env bash

#BSUB -J webpath
#BSUB -o err-out/%J.%I.out
#BSUB -e err-out/%J.%I.err

set -o nounset -o pipefail -o errexit -x

LSB_JOBINDEX=1

source config.sh

trackfile=${WEBFOLDER}${PROJECTDIR}/tracks.txt

if [[ ! -d ${WEBFOLDER}${PROJECTDIR} ]]; then
  mkdir ${WEBFOLDER}${PROJECTDIR}
fi

if [[ -f $trackfile ]]; then
  rm -r $trackfile
  touch $trackfile
fi

NAMEPREFIX=""

for ((i=0; i<${numBarcodes}; i++));
do
	prefix=${BARCODES[$i]}
  NORM_BGFILE=$prefix.norm.bg
  NORM_POSBGFILE=$prefix.norm.pos.bg
  NORM_NEGBGFILE=$prefix.norm.neg.bg
  BIGWIG="$prefix.bw"
  BIGWIGPOS="$prefix.pos.bw"
  BIGWIGNEG="$prefix.neg.bw"
  bedGraphToBigWig $NORM_BGFILE $CHROMSIZES $BIGWIG
  bedGraphToBigWig $NORM_POSBGFILE $CHROMSIZES $BIGWIGPOS
  bedGraphToBigWig $NORM_NEGBGFILE $CHROMSIZES $BIGWIGNEG

	cp $BIGWIG ${WEBFOLDER}${PROJECTDIR}/.
  cp $BIGWIGPOS ${WEBFOLDER}${PROJECTDIR}/.
  cp $BIGWIGNEG ${WEBFOLDER}${PROJECTDIR}/.
	echo 'track type=bigWig name="'${NAMEPREFIX}''${prefix}_all'" description="'${DESCRS[$i]} all'" bigDataUrl="'${SANDBOX}${PROJECTDIR}/${BIGWIG}'" maxHeightPixels=50:35:25 visibility=full color='${COLORS[$i]} >> $trackfile
  echo 'track type=bigWig name="'${NAMEPREFIX}''${prefix}_pos'" description="'${DESCRS[$i]}_pos'" bigDataUrl="'${SANDBOX}${PROJECTDIR}/${BIGWIGPOS}'" maxHeightPixels=50:35:25 visibility=full color='${COLORS[$i]} >> $trackfile
  echo 'track type=bigWig name="'${NAMEPREFIX}''${prefix}_neg'" description="'${DESCRS[$i]}_neg'" bigDataUrl="'${SANDBOX}${PROJECTDIR}/${BIGWIGNEG}'" maxHeightPixels=50:35:25 visibility=full color='${COLORS[$i]} >> $trackfile
done
