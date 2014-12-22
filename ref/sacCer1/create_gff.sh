#Creation of major sacCer1 reference bed

# GFF (saccharomyces_cerevisiae.gff) downloaded on 2014.04.09
# dated 2014.04.08

cd /vol2/home/speach/ref/genomes/sacCer1

tail -n 174966 saccharomyces_cerevisiae.gff > sacCer1.gff
awk '$2 == "SGD" || $2 == "landmark"' sacCer1.gff > sacCer1.full.tab
python gff_parse.py sacCer1.full.tab > sacCer1.full.bed
