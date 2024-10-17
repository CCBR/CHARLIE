import sys
import matplotlib
import numpy
import scipy
#from matplotlib_venn import venn2
#import matplotlib.pyplot as plt

if len(sys.argv)<3:
	print("python %s a_list b_list"%(sys.argv[0]))
	exit()
a_set=set(list(filter(lambda x:x!="",list(map(lambda x:x.strip().split("\t")[0],open(sys.argv[1]).readlines())))))
b_set=set(list(filter(lambda x:x!="",list(map(lambda x:x.strip().split("\t")[0],open(sys.argv[2]).readlines())))))
a_intersect_b=a_set.intersection(b_set)
a_union_b=a_set.union(b_set)
a_only=a_set-b_set
b_only=b_set-a_set
print("Size of a_list=%d"%(len(a_set)))
print("Size of b_list=%d"%(len(b_set)))
print("a interset b=%d"%(len(a_intersect_b)))
print("a union b=%d"%(len(a_union_b)))
print("only a=%d"%(len(a_only)))
print("only b=%d"%(len(b_only)))
if len(sys.argv)==4:
	def write_list_to_file(a_set,filename):
		o=open(filename,'w')
		for g in a_set:
			o.write("%s\n"%(g))
		o.close()
	write_list_to_file(a_intersect_b,"a_intersect_b.lst")
	write_list_to_file(a_union_b,"a_union_b.lst")
	write_list_to_file(a_only,"a_only.lst")
	write_list_to_file(b_only,"b_only.lst")
#venn2(subsets = (len(a_only), len(b_only), len(a_intersect_b)))
#plt.savefig("ab_venn.png")
exit()