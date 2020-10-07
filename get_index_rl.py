import sys
import gzip
from itertools import islice
with gzip.open(sys.argv[1],'r') as fin:
	for line in islice(fin,1,2) :
		r=len(line.strip())

offset=2
rls=[50,75,100,125,150]
b=list(map(lambda x:x-int(r),rls))
c=list(filter(lambda x:x<=(0+offset),b))
print(rls[b.index(max(c))])
