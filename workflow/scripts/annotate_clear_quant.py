# AUTHOR: Vishal N Koparde, Ph.D, CCBR, NCI
# This script annotates quant_txt file from CLEAR workflow with hg38_2_hg19_lookup.txt

import pandas
import re
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
import sys

indexcol=sys.argv[3] # hg38 or mm39

lookupfile=sys.argv[1]
annotations=pandas.read_csv(lookupfile,sep="\t",header=0)
annotations.set_index([indexcol],inplace=True)

quantfile=sys.argv[2]
quant=pandas.read_csv(quantfile,sep="\t",header=None,names=["quant_chrom","quant_start","quant_end","quant_name","quant_score","quant_quant_strand","quant_thickStart","quant_thickEnd","quant_itemRgb","quant_exonCount","quant_exonSizes","quant_exonOffsets","quant_readNumber","quant_circType","quant_geneName","quant_isoformName","quant_index","quant_flankIntron","quant_FPBcirc","quant_FPBlinear","quant_CIRCscore"])
quant[indexcol]=quant.apply(lambda row: row.quant_chrom+":"+str(row.quant_start)+"-"+str(row.quant_end),axis=1)
quant.set_index([indexcol],inplace=True)

x=quant.join(annotations)
x.to_csv(quantfile+'.annotated',sep="\t",header=True)
