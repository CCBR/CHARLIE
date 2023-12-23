import sys
stats=dict()
mreads=int(sys.argv[3]) # minreads
#read junction.filter1
with open(sys.argv[1]) as junction:
    for l in junction.readlines():
        l=l.strip().split("\t")
        if l[0]!=l[3]:continue
        if l[2]!=l[5]:continue
        if l[1]==l[4]:continue
        if int(l[1]) > int(l[4]):
            end = l[1]
            start = l[4]
        else:
            end = l[4]
            start = l[1]
        jid = l[0] + "##" + start + "##" + end
        if not jid in stats:
            stats[jid]=dict()
            stats[jid]["+"]=0
            stats[jid]["-"]=0
        stats[jid][l[2]]+=1
#read back_spliced_junction.filter2.bed
with open(sys.argv[2]) as bsjbed:
    for l in bsjbed.readlines():
        l=l.strip().split("\t")
        if l[1]==l[2]:continue
        jname,count=l[3].split("/")
        if int(count)<mreads:
            continue
        else:
            l[3]=str(count)
        end=int(l[2])+1
        bsjid=l[0] + "##" + l[1] + "##" + str(end)
        if not bsjid in stats:
            strand="."
        else:
            if stats[bsjid]["+"] > stats[bsjid]["-"]:
                strand="+"
            else:
                strand="-"
        l[5]=strand
        print("\t".join(l))