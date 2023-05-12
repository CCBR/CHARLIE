import pysam
import sys
import argparse
import gzip
import os
import time

def get_ctime():
    return time.ctime(time.time())

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
        description="""Extracts PE BSJs from STAR2p output Chimeric BAM file. It also adds
        unique read group IDs to each read. This RID is of the format <chrom>##<start>##<end>
        where the chrom, start and end represent the BSJ the read is depicting.
        ## UPDATE: works for all BAM files ... not just BSJ only 
        """
    )
    #INPUTs
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=str,
        help="Input BAM file")
    parser.add_argument("-s",'--sample_name', dest='samplename', type=str, required=False, default = 'sample1',
        help='Sample Name: SM for RG')
    parser.add_argument('--regions', dest='regions', type=str, required=True,
        help='regions file eg. ref.fa.regions')
    parser.add_argument('--host', dest='host', type=str, required=True,
        help='host name eg.hg38... single value')
    parser.add_argument('--additives', dest='additives', type=str, required=True,
        help='additive name(s) eg.ERCC... comma-separated list... all BSJs in this region are filtered out')
    parser.add_argument('--viruses', dest='viruses', type=str, required=True,
        help='virus name(s) eg.NC_009333.1... comma-separated list')
    parser.add_argument('--prefix', dest='prefix', type=str, required=True,
        help='outfile prefix ... like "linear" or "linear_spliced" etc.')
    #OUTPUTs
    parser.add_argument("--outdir",dest="outdir",required=False,type=str,
        help="Output folder for the individual BAM files.")

    args = parser.parse_args()

    samfile = pysam.AlignmentFile(args.inbam, "rb")
    sequences = list()
    samheader = samfile.header.to_dict()
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
    
    outputbams = dict()
    for h in hosts:
        outbamname = os.path.join(args.outdir,args.samplename+"."+args.prefix+"."+h+".bam")
        outputbams[h] = pysam.AlignmentFile(outbamname, "wb", header = samheader)
    for h in viruses:
        outbamname = os.path.join(args.outdir,args.samplename+"."+args.prefix+"."+h+".bam")
        outputbams[h] = pysam.AlignmentFile(outbamname, "wb", header = samheader)
    
    for read in samfile.fetch():
        chrom=read.reference_name
        regionname=seqname2regionname[chrom]
        if regionname in hosts or regionname in viruses:
            outputbams[regionname].write(read)
    samfile.close()
    for o in outputbams.values():
        o.close()



if __name__ == "__main__":
    main()