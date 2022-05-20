#!/usr/bin/env python3
import sys
import textwrap

usage_txt=textwrap.dedent("""\
	Description:
	The script collapses bed entries, ie, if the bed file has repeated 
	regions but with different names, then they are all collaped into
	a single bed entry and the names are reported as a comma separated
	list in the 4th column
	Usage:
	python3 {} <inputBed6> <outputBed6>
	@Parameters:
	1. <inputBed6>: BED6 file that needs to be collapsed by name
	2. <outputBed6>: BED6 collaped output file
""".format(__file__))


if len(sys.argv)!=3:
	exit(usage_txt)

with open(sys.argv[1]) as f:
	inputBedLines=f.readlines()

names=dict()
for l in inputBedLines:
	l=l.strip().split("\t")
	tmp=[l[0],l[1],l[2],l[5]]
	region_id="##".join(tmp)
	if not region_id in names:
		names[region_id]=list()
	names[region_id].append(l[3])

outbed = open(sys.argv[2],'w')
for region_id,name in names.items():
	tmp=region_id.split("##")
	namelist=",".join(name)
	tmp.insert(3,namelist)
	tmp.insert(4,"0")
	outbed.write("\t".join(tmp)+"\n")
outbed.close()
