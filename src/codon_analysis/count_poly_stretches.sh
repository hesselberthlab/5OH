#! /usr/bin/env bash

aminos=(DDDDD EEEEE)
ref=~/ref/genomes/sacCer1/

for amino in $aminos; do
    echo $amino
    bioawk -c fastx '{print $name,'\t',$seq}' ${ref}orf_trans.fasta \
        | grep $amino \
        | wc | echo
done
    

