table mRNAs 
"custom mRNAs for sacCer1"
(
string  chrom;      "Reference sequence chromosome or scaffold"
uint    chromStart; "Start position of feature on chromosome"
uint    chromEnd;   "End position of feature on chromosome"
string  name;       "Name of gene"
uint    score;      "Score"
char[1] strand;     "+ or - for strand"
uint    thickStart; "Coding region start"
uint    thickEnd;   "Coding region end"
uint reserved;     "Used as itemRgb as of 2004-11-22"
int blockCount;    "Number of blocks"
int[blockCount] blockSizes; "Comma separated list of block sizes"
int[blockCount] chromStarts; "Start positions relative to chromStart"
string  sysname;    "systematic name"
)
