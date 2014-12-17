genes=(SPT6 FUN12 CBF5 RPN2 ADE8 ADK1 ADH1 RPS31)

echo -e "gene\tstrand\tlength\tcoords"

for gene in ${genes[@]}; do
    grep $gene ~/ref/genomes/sacCer1/sacCer1.mRNA.bed \
    | awk 'BEGIN{OFS="\t"} {print $4, $6, $3-$2, $1 ":" $2 "-" $3}'
done
