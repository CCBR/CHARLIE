import sys
collection=dict()
for f in sys.stdin:
        f=f.strip().split("\t")
        circid="##".join([f[0],f[1],f[2],f[5]])
        if not circid in collection:
                collection[circid]=dict()
                collection[circid]['fullline']=f
                collection[circid]['count']=int(f[4])
        else:
                collection[circid]['count']+=int(f[4])
#header=["chrom","start","end","name","n_reads","strand","n_uniq","uniq_bridges","best_qual_left","best_qual_right","tissues","tiss_counts","edits","anchor_overlap","breakpoints","signal","strandmatch","category"]
#print("\t".join(header))
count=0
for k,v in collection.items():
        count+=1
        x=v['fullline']
        x[3]=str(count)
        x[4]=str(v['count'])
        print("\t".join(x))