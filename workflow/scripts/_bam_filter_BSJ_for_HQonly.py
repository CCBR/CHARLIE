import argparse
import pandas as pd
import pysam
import os

def read_regions(regionsfile,host,additives,viruses):
    host=host.split(",")
    additives=additives.split(",")
    viruses=viruses.split(",")
    infile=open(regionsfile,'r')
    regions=dict()
    for l in infile.readlines():
        l = l.strip().split("\t")
        region_name=l[0]
        regions[region_name]=dict()
        regions[region_name]['sequences']=dict()
        if region_name in host:
            regions[region_name]['host_additive_virus']="host"
        elif region_name in additives:
            regions[region_name]['host_additive_virus']="additive"
        elif region_name in viruses:
            regions[region_name]['host_additive_virus']="virus"
        else:
            exit("%s has unknown region. Its not a host or a additive or a virus!!")
        sequence_names=l[1].split()
        for s in sequence_names:
            regions[region_name]['sequences'][s]=1
    return regions

def _get_host_additive_virus(regions,seqname):
    for k,v in regions.items():
        if seqname in v['sequences']:
            return v['host_additive_virus']
    else:
        exit("Sequence: %s does not have a region."%(seqname))

def _get_regionname_from_seqname(regions,seqname):
    for k,v in regions.items():
        if seqname in v['sequences']:
            return k
    else:
        exit("Sequence: %s does not have a region."%(seqname))     

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
        description="""RG is created from info from the counts table.
        This RG is used to extract reads from inbam and save them.
        """
    )
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=str,
        help="BSJ bam with RG set")
    parser.add_argument('-t','--sample_counts_table', dest='countstable', type=str, required=True,
        help='final all sample counts matrix')	# get coordinates of the circRNA
    # parser.add_argument("-o","--outbam",dest="outbam",required=True,type=argparse.FileType('w'),
    #     help="Output bam file ... both strands")
    parser.add_argument("-s",'--sample_name', dest='samplename', type=str, required=False, default = 'sample1',
        help='Sample Name')
    parser.add_argument("-o","--outbam",dest="outbam",required=True,type=str,
        help="Output bam file ... both strands")
    parser.add_argument('--regions', dest='regions', type=str, required=True,
        help='regions file eg. ref.fa.regions')
    parser.add_argument('--host', dest='host', type=str, required=True,
        help='host name eg.hg38... single value')
    parser.add_argument('--additives', dest='additives', type=str, required=True,
        help='additive name(s) eg.ERCC... comma-separated list... all BSJs in this region are filtered out')
    parser.add_argument('--viruses', dest='viruses', type=str, required=True,
        help='virus name(s) eg.NC_009333.1... comma-separated list')
    args = parser.parse_args()		

    indf = pd.read_csv(args.countstable,sep="\t",header=0,compression='gzip')
    indf = indf.loc[indf['HQ']=="Y"]

    RGlist = dict()
    for index,row in indf.iterrows():
        jid = row['chrom']+"##"+str(row['start'])+"##"+str(row['end'])
        RGlist[jid]=1
    print("Number of RGs: ",len(RGlist))

    samfile = pysam.AlignmentFile(args.inbam, "rb")
    samheader = samfile.header.to_dict()

    sequences = list()
    for v in samheader['SQ']:
        sequences.append(v['SN'])
    seqname2regionname=dict()
    hosts=set()
    viruses=set()
    regions = read_regions(regionsfile=args.regions,host=args.host,additives=args.additives,viruses=args.viruses)
    for s in sequences:
        hav = _get_host_additive_virus(regions,s)
        if hav == "host":
            hostname = _get_regionname_from_seqname(regions,s)
            seqname2regionname[s]=hostname
            hosts.add(hostname)
        if hav == "virus":
            virusname = _get_regionname_from_seqname(regions,s)
            seqname2regionname[s]=virusname
            viruses.add(virusname)


    outbam = pysam.AlignmentFile(args.outbam, "wb", template=samfile)
    outputbams = dict()
    outdir = os.path.dirname(args.outbam)
    for h in hosts:
        outbamname = os.path.join(outdir,args.samplename+"."+h+".HQ_only.BSJ.bam")
        outputbams[h] = pysam.AlignmentFile(outbamname, "wb", header = samheader)

    for v in viruses:
        outbamname = os.path.join(outdir,args.samplename+"."+v+".HQ_only.BSJ.bam")
        outputbams[v] = pysam.AlignmentFile(outbamname, "wb", header = samheader)


    for read in samfile.fetch():
        rg = read.get_tag("RG")
        rg = rg.split("##")
        rg = rg[:len(rg)-1]
        rg = "##".join(rg)
        if rg in RGlist:
            regionname=_get_regionname_from_seqname(regions,read.reference_name)
            if regionname in hosts:
                outputbams[regionname].write(read)
            if regionname in viruses:
                outputbams[regionname].write(read)
            outbam.write(read)
    samfile.close()
    outbam.close()
    for k,v in outputbams.items():
        v.close()


if __name__ == "__main__":
    main()
