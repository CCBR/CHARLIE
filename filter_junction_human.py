import sys
for i in open(sys.argv[1]).readlines():
	j=i.split("\t")
	if (j[0]!="chrKSHV") and (j[3]==j[0]) :
		print(i.strip())

