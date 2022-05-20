# AUTHOR: Vishal N Koparde, Ph.D, CCBR, NCI
# This script annotates quant.txt file from CLEAR workflow with hg38_2_hg19_lookup.txt

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
quant=pandas.read_csv(quantfile,sep="\t",header=None,names=["chrom","start","end","name","score","quant_strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes","exonOffsets","readNumber","circType","geneName","isoformName","index","flankIntron","FPBcirc","FPBlinear","CIRCscore"])
quant[indexcol]=quant.apply(lambda row: row.chrom+":"+str(row.start)+"-"+str(row.end),axis=1)
quant.set_index([indexcol],inplace=True)

x=quant.join(annotations)
x.to_csv(quantfile+'.annotated',sep="\t",header=True)
