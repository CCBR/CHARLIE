import argparse
import pandas

debug=0

def get_attributes(attstr):
	att = dict()
	attlist = attstr.strip().split(";")
	if debug==1: print(attstr)
	if debug==1: print(attlist)
	for item in attlist:
		x = item.strip()
		if debug==1: print(x)
		x = x.replace("\"","")
		if debug==1: print(x)
		x = x.split()
		if debug==1: print(x)
		if len(x)!=2: continue
		key = x.pop(0)
		key = key.replace(":","")
		value = " ".join(x)
		value = value.replace(":","_")
		att[key] = value
	return att

def get_attstr(att):
	strlist=[]
	for k,v in att.items():
		s = "%s \"%s\""%(k,v)
		strlist.append(s)
	attstr = "; ".join(strlist)
	return attstr+";"

parser = argparse.ArgumentParser(description='fix gtf file')
parser.add_argument('--ingtf', dest='ingtf', type=str, required=True,
		                    help='input gtf file')
parser.add_argument('--outgtf', dest='outgtf', type=str, required=True,
		                    help='output gtf file')
args = parser.parse_args()

gene_id_2_gene_name = dict()

with open(args.ingtf, 'r') as ingtf:
	for line in ingtf:
		if line.startswith("#"): continue
		line = line.strip()
		line = line.split("\t")
		if len(line) != 9:
			print(line)
			exit("ERROR ... line does not have 9 items!")
		attributes = get_attributes(line[8])
		if debug==1: print(line)
		if debug==1: print(attributes)
		if not attributes["gene_id"] in gene_id_2_gene_name:
			if "gene_name" in attributes:
				gene_id_2_gene_name[attributes["gene_id"]] = attributes["gene_name"]
			else:
				gene_id_2_gene_name["gene_id"] = attributes["gene_id"]

with open("gene_id_2_gene_name.tsv",'w') as tmp:
	for k,v in gene_id_2_gene_name.items():
		tmp.write("%s\t%s\n"%(k,v))

with open(args.ingtf,'r') as ingtf, open(args.outgtf,'w') as outgtf:
	for line in ingtf:
		if line.startswith("#"): 
			outgtf.write(line)
			continue
		line = line.strip()
		line = line.split("\t")
		attributes = get_attributes(line[8])
		if not "gene_name" in attributes:
			if not "gene_id" in attributes:
				print(line)
				print(attributes)
				exit("ERROR in this line!")
			if not attributes["gene_id"] in gene_id_2_gene_name:
				print(line)
				print(attributes)
				print(attributes["gene_id"])
				exit("ERROR2 in this line!")
			attributes["gene_name"] = gene_id_2_gene_name[attributes["gene_id"]]
		line[8] = get_attstr(attributes)
		outgtf.write("\t".join(line)+"\n")

		
