import sys
import argparse
import os

def my_bool(s):
    return s != 'False'

parser = argparse.ArgumentParser(description='apply junction filters, input stdin and output stdout')
parser.add_argument('--regions', dest='regions', type=str, required=True, metavar="absolute path to regions file",
                    help='regions file')
parser.add_argument('--filter1regions', dest='filter1regions', type=str, required=True, metavar="eg. \"hg38,ERCC,rRNA\"",
                    help='comma separated list of regions to apply filter1 on ... filter2 is applied to all other regions')
parser.add_argument('--filter1_noncanonical', dest='filter1_noncanonical', default=True, type=my_bool, required=True, metavar="\"True/False\"",
                    help='apply canonical filter on filter1')
parser.add_argument('--filter1_unannotated', dest='filter1_unannotated', default=True, type=my_bool, required=True, metavar="\"True/False\"",
                    help='apply unannotated filter on filter1')
parser.add_argument('--filter2_noncanonical', dest='filter2_noncanonical', default=False, type=my_bool, required=True, metavar="\"True/False\"",
                    help='apply canonical filter on filter2')
parser.add_argument('--filter2_unannotated', dest='filter2_unannotated', default=False, type=my_bool, required=True, metavar="\"True/False\"",
                    help='apply unannotated filter on filter2')					
args = parser.parse_args()

chr2region=dict()
regions=list()
x = open(args.regions)
for r in x.readlines():
	r = r.strip().split("\t")
	regions.append(r[0])
	for c in r[1].split():
		chr2region[c]=r[0]
x.close()

region2filter=dict()
for x in regions:
	region2filter[x]=2 # apply filter2 to everything

filter1regions=args.filter1regions
for f in filter1regions.split(","):
	f = f.strip()
	if not f in region2filter:
		exit("Region "+f+" not defined!")
	region2filter[f]=1 # change filter from filter2 to filter1

# cat {input} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort -k1,1 -k2,2n|uniq > {output.pass1sjtab}
for line in sys.stdin:
	l=line.split("\t")
	f=region2filter[chr2region[l[0]]]
	if f==1:
		if args.filter1_noncanonical:
			if not int(l[4])>0:
				continue
		if args.filter1_unannotated:
			if not int(l[5])==1:
				continue
	elif f==2:
		if args.filter2_noncanonical:
			if not int(l[4])>0:
				continue
		if args.filter2_unannotated:
			if not int(l[5])==1:
				continue
	sys.stdout.write(line)
	# exit()

