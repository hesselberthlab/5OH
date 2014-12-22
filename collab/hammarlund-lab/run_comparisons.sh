#! /usr/bin/env bash

<<DOC
Compare changepoint combinations in 4 Hammarlund conditions
DOC

PIPELINE="/vol2/home/speach/devel/5OH/src/collab/hammarlund-lab/"
CONDITIONS=(CNTRL_0 CNTRL_TM RTCB_0 RTCB_TM)

len=${#CONDITIONS[@]}

for (( i=0; i<${len}; i++));
do
    for (( j=0; j<${len}; j++));
    do
        if [[ $i == $j ]]; then
            continue
        fi

        map_a=${CONDITIONS[$i]}*exonmappings.tab
        map_b=${CONDITIONS[$j]}*exonmappings.tab
        name_a=$(basename ${CONDITIONS[$i]} _0)
        name_b=$(basename ${CONDITIONS[$j]} _0)
        compname=${name_a}.vs.${name_b}
        
        bash ${PIPELINE}compare_changepoints_map.sh $map_a $map_b $compname
    done
done