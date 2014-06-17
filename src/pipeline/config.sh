PROJECTID=5OH-methods

# this sets the RESULT directory to "common-debug"
DEBUG="-debug"

# XXX SAMPLES, DESCRIPS, and COLORS must be in same order

# Rep1
SAMPLES=(SP5
          SP6
          SP7
          SP8
          SP9
          SP10
          SP11
          SP12)

NUM_SAMPLES=${#SAMPLES[@]}

DESCRIPS=("del-XRN1, +Tm - Rep1"
          "del-XRN1, DMSO - Rep1"
          "WT, +Tm - Rep1"
          "WT, DMSO - Rep1"
          "del-XRN1, +Tm, +SAP - Rep1"
          "del-XRN1, DMSO, +SAP - Rep1"
          "WT, +Tm, +SAP - Rep1"
          "WT, DMSO, +SAP - Rep1")

RED="215,25,28"
ORANGE="253,174,97"
GREEN="0,100,0"
BLUE="43,131,186"

COLORS=($RED
        $ORANGE
        $GREEN
        $BLUE
        $RED
        $ORANGE
        $GREEN
        $BLUE)

DATA=$HOME/projects/5OH/data/Rep1
BIN=$HOME/devel/5OH
RSCRIPTS=$HOME/devel/5OH/src/R
UMI=NNNNNNNN

# Variable alignments created by using different arguments to bowtie
# For yeast, two copies of rDNA, so need at least 2 alignments to capture
ALIGN_MODES=("uniq" "two" "all")
ALIGN_ARGS=("-m 1" "-m 2" "--all")