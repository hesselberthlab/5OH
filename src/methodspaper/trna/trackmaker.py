'''
Create UCSC genome browser tracks which include sample descriptions
'''

import os
import sys

results_dir = sys.argv[1]

COLOR_DICT = {"red":"215,25,28", "green":"0,100,0",
          "blue":"43,131,186", "orange":"253,174,97"}

# Manually enter experiment-specific sample information
# Dictionary keyed by sample ID with values as a tuple of: (sample
# description, color for track)
SAMPLE_DATA = {"SP5":("delXrn1_Tm_Rep1","red"),
               "SP6":("delXrn1_DMSO_Rep1","orange"),
               "SP7":("WT_Tm_Rep1", "green"),
               "SP8":("WT_DMSO_Rep1", "blue"),
               "SP9":("delXrn1_Tm_+SAP_Rep1","red"),
               "SP10":("delXrn1_DMSO_+SAP_Rep1","orange"),
               "SP11":("WT_Tm_+SAP_Rep1", "green"),
               "SP12":("WT_DMSO_+SAP_Rep1", "blue"),
               "SP13":("delXrn1_Tm_Rep2","red"),
               "SP27":("delXrn1_DMSO_Rep2","orange"),
               "SP14":("WT_Tm_Rep2", "green"),
               "SP15":("WT_DMSO_Rep2", "blue"),
               "SP16":("delXrn1_Tm_+SAP_Rep2","red"),
               "SP17":("delXrn1_DMSO_+SAP_Rep2","orange"),
               "SP18":("WT_Tm_+SAP_Rep2", "green"),
               "SP19":("WT_DMSO_+SAP_Rep2", "blue"),
               "SP28":("WT_RiboZero_Rep1", "red"),
               "SP20":("WT_Pre-Frag_Rep1", "orange"),
               "SP29":("WT_Pre-Frag_Rep2", "orange")}

#nah, just a few
SAMPLE_DATA ={"SP8":("WT_DMSO_Rep1_smallshit", "blue")}

SANDBOX = "http://amc-sandbox.ucdenver.edu/~speach/"
WEBURL = "%s%s" % (SANDBOX, results_dir)
PUBLICHTMLDIR = "/vol2/home/speach/public_html/%s" % results_dir

# Cycle through all bigwig files in directory and print tracks 
for filename in os.listdir(PUBLICHTMLDIR):
    filefields = filename.strip().split(".")

    # select bigwigs of CPMs (change second argument to "counts"
    # to get raw counts) (change third argument for strand stuff)
    if filefields[-1] == "bw" and filefields[-2] == "CPMs":
        sample, assembly, alignment = filefields[:3]
        assembly = assembly.split("-")[1]
        alignment = alignment.split("-")[1]

        sample_descr = SAMPLE_DATA[sample][0]
        sample_color = SAMPLE_DATA[sample][1]
        sample_color_hex = COLOR_DICT[sample_color]

        track = "track type=bigWig"
        name = 'name="%s_%s-%s"' % (sample, assembly, alignment)
        descr = 'description="%s"' % sample_descr
        url = 'bigDataUrl="%s%s"' % (WEBURL, filename)
        height = "maxHeightPixels=50:35:25"
        vis = "visibility=full"
        color = "color=%s" % sample_color_hex

        if alignment == "small" and filefields[-3] == "all":
            print " ".join([track, name, descr, url, height, vis, color])
