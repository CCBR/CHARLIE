import sys
def get_id(s,whatid):
	s=s.split()
	for i,j in enumerate(s):
		if j==whatid:
			r=s[i+1]
	r=r.replace('"','')
	r=r.replace(';','')
	return r

		

gtffile=sys.argv[1]
transcript2gene=dict()
for i in open(gtffile).readlines():
	if i.startswith("#"):
		continue
	i=i.strip().split("\t")
	if i[2]!="transcript":
		continue
	gid=get_id(i[8],"gene_id")
	tid=get_id(i[8],"transcript_id")
#	print("%s\t%s"%(tid,gid))
	transcript2gene[tid]=gid

for i in open(sys.argv[2]).readlines():
	j=i.strip().split("\t")
	x=[]
	tid=j.pop(0)
	gid=transcript2gene[tid]
	x.append(gid)
	x.append(tid)
	x.extend(j)
	print("\t".join(x))
