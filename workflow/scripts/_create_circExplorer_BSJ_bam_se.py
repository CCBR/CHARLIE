import pysam
import sys
import argparse
import gzip
import os
import time

def get_ctime():
    return time.ctime(time.time())

"""

This script first validates each read to be "valid" BSJ read and then splits a BSJ bam file by strand into:
1. Plus strand only bam file
2. Minus strand only bam file
3. BSJ bed file with score(number of reads supporting the BSJ) and strand information
Logic (for SE reads):
Each BSJ is represented by a 2 alignments in the output BAM file.
Alignments 1 and 2 are split alignment of read1 at two distinct loci on the same reference 
chromosome.
These alignments are grouped together by the "HI" tags in SAM file. For example, all 2 
alignments for the same BSJ will have the same "HI" value... something like "HI:i:1".
BSJ alignment sam bitflag combinations can have 4 different possibilities, 2 from sense strand
and 2 from anti-sense strand:
1. 0, 2048 or 2048, 0
2. 256, 2304 or 2304, 256
3. 16, 2064 or 2064, 16
4. 256, 2320 or 2320, 256
"""


class BSJ:
    def __init__(self):
        self.chrom=""
        self.start=""
        self.end=""
        self.score=0
        self.name="."
        self.strand="U"
        self.bitids=set()
        self.rids=set()
        
    def plusone(self):
        self.score+=1
    
    def set_strand(self,strand):
        self.strand=strand
    
    def set_chrom(self,chrom):
        self.chrom=chrom

    def set_start(self,start):
        self.start=start

    def set_end(self,end):
        self.end=end
    
    def append_bitid(self,bitid):
        self.bitids.add(bitid)

    def append_rid(self,rid):
        self.rids.add(rid)
        
    def write_out_BSJ(self,outbed):
        t=[]
        t.append(self.chrom)
        t.append(str(self.start))
        t.append(str(self.end))
        t.append(self.name)
        t.append(str(self.score))
        t.append(self.strand)
        t.append(",".join(self.bitids))
        t.append(",".join(self.rids))
        outbed.write("\t".join(t)+"\n")		

    def update_score_and_found_count(self,junctions_found):
        self.score = len(self.rids)
        jid = self.chrom + "##" + str(self.start) + "##" + str(int(self.end)-1) + "##" + self.strand
        junctions_found[jid]+=self.score
        
class Readinfo:
    def __init__(self,readid,rname):
        self.readid=readid
        self.refname=rname
        # self.alignments=list()
        self.bitflags=list()
        self.bitid=""
        self.strand="."
        self.start=-1
        self.end=-1
        self.refcoordinates=dict()
        self.isread1=dict()
        self.isreverse=dict()
        self.issecondary=dict()
        self.cigarstrs=dict()
        self.issupplementary=dict()
    
    def __str__(self):
        s = "readid: %s"%(self.readid)
        s = "%s\tbitflags: %s"%(s,self.bitflags)
        s = "%s\tisreverse: %s"%(s,self.isreverse)
        s = "%s\tbitid: %s"%(s,self.bitid)
        return s

    def set_refcoordinates(self,bitflag,refpos):
        self.refcoordinates[bitflag]=refpos
    
    def set_cigarstr(self,bitflag,cigarstr):
        self.cigarstrs[bitflag]=cigarstr
    
    def set_read1_reverse_secondary_supplementary(self,bitflag,read):
        if read.is_read1:
            self.isread1[bitflag]="Y"
        else:
            self.isread1[bitflag]="N"
        if read.is_reverse:
            self.isreverse[bitflag]="Y"
        else:
            self.isreverse[bitflag]="N"
        if read.is_secondary:
            self.issecondary[bitflag]="Y"
        else:
            self.issecondary[bitflag]="N"
        if read.is_supplementary:
            self.issupplementary[bitflag]="Y"
        else:
            self.issupplementary[bitflag]="N"
    
    # def append_alignment(self,read):
    #     self.alignments.append(read)
    
    def append_bitflag(self,bf):
        self.bitflags.append(bf)
    
    # def extend_ref_positions(self,refcoords):
    # 	self.refcoordinates.extend(refcoords)
    
    def generate_bitid(self):
        bitlist=sorted(self.bitflags)
        self.bitid="##".join(list(map(lambda x:str(x),bitlist)))
# 		self.bitid=str(bitlist[0])+"##"+str(bitlist[1])+"##"+str(bitlist[2])
    
    def get_strand(self):
        if self.bitid=="0##2048":
            self.strand="-"
        elif self.bitid=="256##2304":
            self.strand="-"
        elif self.bitid=="16##2064":
            self.strand="+"
        elif self.bitid=="272##2320":
            self.strand="+"		
        else:
            self.strand="U"

    def validate_BSJ_read(self,junctions):
        """
        Checks if read is truly a BSJ originitor.
        """
        if len(self.bitid.split("##"))==2:
            if not self.bitid in ["0##2048", "16##2064", "256##2304", "272##2320"]:
                return False
            count=0
            refcoords=self.refcoordinates
            for k,v in refcoords.items():
                count+=1
                refcoords[k]=sorted(v)
                if count==1:
                    astart=refcoords[k][0]
                    aend=refcoords[k][-1]
                if count==2:
                    bstart=refcoords[k][0]
                    bend=refcoords[k][-1]
            chrom = self.refname
            possiblejid=chrom+"##"+str(astart)+"##"+str(bend)+"##"+self.strand
            possiblejid2=chrom+"##"+str(bstart)+"##"+str(aend)+"##"+self.strand
            # exit()
            if possiblejid in junctions:
                self.start = astart
                self.end = str(int(bend) + 1)    # this will be added to the BED file
                return True
            if possiblejid2 in junctions:
                self.start = bstart
                self.end = str(int(aend) + 1)    # this will be added to the BED file
                return True            
        else:
            return False
    
    def get_bsjid(self):
        t=[]
        t.append(self.refname)
        t.append(str(self.start))
        t.append(str(self.end))
        t.append(self.strand)
        return "##".join(t)
    
    # def write_out_reads(self,outbam):
    #     for r in self.alignments:
    #         outbam.write(r)
        
            
def get_uniq_readid(r):
    rname=r.query_name
    hi=r.get_tag("HI")
    rid=rname+"##"+str(hi)
    return rid

def get_bitflag(r):
    bitflag=str(r).split("\t")[1]
    return int(bitflag)

def _bsjid2chrom(bsjid):
    x=bsjid.split("##")
    return x[0]

def _bsjid2jid(bsjid):
    x=bsjid.split("##")
    chrom=x[0]
    start=x[1]
    end=str(int(x[2])-1)
    jid="##".join([chrom,start,end])
    return jid,chrom

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
        description="""Extracts SE BSJs from STAR2p output Chimeric BAM file. It also adds
        unique read group IDs to each read. This RID is of the format <chrom>##<start>##<end>
        where the chrom, start and end represent the BSJ the read is depicting.
        """
    )
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=str,
        help="Input Chimeric-only STAR2p BAM file")
    parser.add_argument("-s",'--sample_name', dest='samplename', type=str, required=False, default = 'sample1',
        help='Sample Name: SM for RG')
    parser.add_argument("-l",'--library', dest='library', type=str, required=False, default = 'lib1',
        help='Sample Name: LB for RG')
    parser.add_argument("-f",'--platform', dest='platform', type=str, required=False, default = 'illumina',
        help='Sample Name: PL for RG')
    parser.add_argument("-u",'--unit', dest='unit', type=str, required=False, default = 'unit1',
        help='Sample Name: PU for RG')
    parser.add_argument('-t','--sample_counts_table', dest='countstable', type=str, required=True,
        help='circExplore per-sample counts table')	# get coordinates of the circRNA
    parser.add_argument("-p","--plusbam",dest="plusbam",required=True,type=argparse.FileType('w'),
        help="Output plus strand bam file")
    parser.add_argument("-m","--minusbam",dest="minusbam",required=True,type=argparse.FileType('w'),
        help="Output plus strand bam file")
    parser.add_argument("-o","--outbam",dest="outbam",required=True,type=argparse.FileType('w'),
        help="Output bam file ... both strands")
    parser.add_argument("--outputhostbams",dest="outputhostbams",required=False,action='store_true', default=False,
        help="Output individual host BAM files")
    parser.add_argument("--outputvirusbams",dest="outputvirusbams",required=False,action='store_true', default=False,
        help="Output individual virus BAM files")
    parser.add_argument("--outdir",dest="outdir",required=False,type=str,
        help="Output folder for the individual BAM files (required only if --outputhostbams or --outputvirusbams is used).")
    parser.add_argument("-b","--bed",dest="bed",required=True,type=str,
        help="Output BSJ bed.gz file (with strand info)")
    parser.add_argument("-j","--junctionsfound",dest="junctionsfound",required=True,type=argparse.FileType('w', encoding='UTF-8'),
        help="Output TSV file with counts of junctions expected vs found")
    parser.add_argument('--regions', dest='regions', type=str, required=True,
        help='regions file eg. ref.fa.regions')
    parser.add_argument('--host', dest='host', type=str, required=True,
        help='host name eg.hg38... single value')
    parser.add_argument('--additives', dest='additives', type=str, required=True,
        help='additive name(s) eg.ERCC... comma-separated list... all BSJs in this region are filtered out')
    parser.add_argument('--viruses', dest='viruses', type=str, required=True,
        help='virus name(s) eg.NC_009333.1... comma-separated list')
    args = parser.parse_args()		
    samfile = pysam.AlignmentFile(args.inbam, "rb")
    samheader = samfile.header.to_dict()
    samheader['RG']=list()
# 	bsjfile = open(args.bed,"w")
    junctionsfile = open(args.countstable,'r')
    junctions=dict()
    junctions_found=dict()
    print("%s | Reading...junctions!..."%(get_ctime()))
    for l in junctionsfile.readlines():
        if "read_count" in l: continue
        l = l.strip().split("\t")
        chrom = l[0]
        start = l[1]
        end = str(int(l[2])-1)
        strand = l[3]
        jid = chrom+"##"+start+"##"+end+"##"+strand                     # create a unique junction ID for each line in the BSJ junction file and make it the dict key ... easy for searching!
        samheader['RG'].append({'ID':jid, 'LB':args.library, 'PL':args.platform, 'PU':args.unit,'SM':args.samplename})
        junctions[jid] = int(l[4])
        junctions_found[jid] = 0
    junctionsfile.close()
    # print(junctions)
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
    print("%s | Done reading %d junctions."%(get_ctime(),len(junctions)))


    bigdict=dict()
    print("%s | Reading...alignments!..."%(get_ctime()))
    count=0
    count2=0
    for read in samfile.fetch():
        count+=1
        satag=read.get_tag("SA")
        satagchrids=list(map(lambda x:samfile.get_tid(x),list(filter(lambda x:x!='',list(map(lambda x:x.split(",")[0],satag.split(";")))))))
        if not read.reference_id in satagchrids: continue    # specific for SE as read.next_reference_id is -1 for SE
        count2+=1
        rid=get_uniq_readid(read)                           # add the HI number to the readid
        if debug:print(rid)
        if not rid in bigdict:
            bigdict[rid]=Readinfo(rid,read.reference_name)
        # bigdict[rid].append_alignment(read)                 # since rid has HI number included ... this separates alignment by HI
        bitflag=get_bitflag(read)
        if debug:print(bitflag)
        bigdict[rid].append_bitflag(bitflag)                # each rid can have upto 3 lines in the BAM with each having its own bitflag ... collect all bitflags in a list here 
        refpos=list(filter(lambda x:x!=None,read.get_reference_positions(full_length=True)))
        # if debug:print(refpos)
        bigdict[rid].set_refcoordinates(bitflag,refpos)     # maintain a list of reference coordinated that are "aligned" for each bitflag in each rid alignment
        bigdict[rid].set_cigarstr(bitflag,read.cigarstring)
        bigdict[rid].set_read1_reverse_secondary_supplementary(bitflag,read)
        if debug:print(bigdict[rid])
    print("%s | Done reading %d chimeric alignments. [%d same chrom chimeras]"%(get_ctime(),count,count2))
    if debug:
        for rid in bigdict.keys():
            print(">>>%s\t%s\t%s\t%s"%(rid,bigdict[rid].isreverse,bigdict[rid].cigarstrs,bigdict[rid].refcoordinates))
    samfile.reset()

    print("%s | Writing BAMs"%(get_ctime()))
    plusfile = pysam.AlignmentFile(args.plusbam, "wb", header = samheader)
    minusfile = pysam.AlignmentFile(args.minusbam, "wb", header = samheader)
    outfile = pysam.AlignmentFile(args.outbam, "wb", header = samheader)
    outputbams = dict()
    if args.outputhostbams:
        for h in hosts:
            outbamname = os.path.join(args.outdir,args.samplename+"."+h+".BSJ.bam")
            outputbams[h] = pysam.AlignmentFile(outbamname, "wb", header = samheader)
    if args.outputvirusbams:
        for v in viruses:
            outbamname = os.path.join(args.outdir,args.samplename+"."+v+".BSJ.bam")
            outputbams[v] = pysam.AlignmentFile(outbamname, "wb", header = samheader)     
    bsjdict=dict()
    bitid_counts=dict()
    lenoutputbams = len(outputbams)
    for read in samfile.fetch():
        satag=read.get_tag("SA")
        satagchrids=list(map(lambda x:samfile.get_tid(x),list(filter(lambda x:x!='',list(map(lambda x:x.split(",")[0],satag.split(";")))))))
        if not read.reference_id in satagchrids: continue    # specific for SE as read.next_reference_id is -1 for SE
        rid=get_uniq_readid(read)
        if rid in bigdict:    
            bigdict[rid].generate_bitid()                       # separate all bitflags for the same rid with ## and create a unique single bitflag ... bitflags are pre-sorted
            if debug:print(bigdict[rid])                        
            bigdict[rid].get_strand()                           # use the unique aggregated bitid to extract the strand information ... all possible cases are explicitly covered
            if not bigdict[rid].validate_BSJ_read(junctions=junctions): # ensure that the read alignments leftmost and rightmost coordinates match with one of the BSJ junctions... if yes then that rid represents a BSJ. Also add start and end to the BSJ object
                continue
            # bigdict[rid].get_start_end()
            # print(bigdict[rid])
            bsjid=bigdict[rid].get_bsjid()
            chrom=_bsjid2chrom(bsjid)
            # jid,chrom=_bsjid2jid(bsjid)
            read.set_tag("RG", bsjid, value_type="Z")
            if bigdict[rid].strand=="+":
                plusfile.write(read)
            if bigdict[rid].strand=="-":
                minusfile.write(read)
            outfile.write(read)
            if lenoutputbams != 0:
                regionname=_get_regionname_from_seqname(regions,chrom)
                if regionname in hosts and args.outputhostbams:
                    outputbams[regionname].write(read)
                if regionname in viruses and args.outputvirusbams:
                    outputbams[regionname].write(read)
            if not bsjid in bsjdict:
                bsjdict[bsjid]=BSJ()
                bsjdict[bsjid].set_chrom(bigdict[rid].refname)
                bsjdict[bsjid].set_start(bigdict[rid].start)
                bsjdict[bsjid].set_end(bigdict[rid].end)
                bsjdict[bsjid].set_strand(bigdict[rid].strand)
            # bsjdict[bsjid].plusone()
            bsjdict[bsjid].append_bitid(bigdict[rid].bitid)
            if not bigdict[rid].bitid in bitid_counts:
                bitid_counts[bigdict[rid].bitid]=0
            bitid_counts[bigdict[rid].bitid]+=1
            bsjdict[bsjid].append_rid(rid)
    plusfile.close()
    minusfile.close()
    samfile.close()
    outfile.close()
    if lenoutputbams != 0:
        for k,v in outputbams.items():
            v.close()
    print("%s | Done!"%(get_ctime()))	
    for b in bitid_counts.keys():
        print(b,bitid_counts[b])
    print("%s | Writing BED"%(get_ctime()))
    with gzip.open(args.bed,'wt') as bsjfile:
        for bsjid in bsjdict.keys():
            bsjdict[bsjid].update_score_and_found_count(junctions_found)
            bsjdict[bsjid].write_out_BSJ(bsjfile)
    bsjfile.close()

    args.junctionsfound.write("#chrom\tstart\tend\tstrand\texpected_BSJ_reads\tfound_BSJ_reads\n")
    for jid in junctions.keys():
        x=jid.split("##")
        chrom=x[0]
        start=int(x[1])
        end=int(x[2])+1
        strand=x[3]
        args.junctionsfound.write("%s\t%d\t%d\t%s\t%d\t%d\n"%(chrom,start,end,strand,junctions[jid],junctions_found[jid]))
    args.junctionsfound.close()
    print("%s | ALL Done!"%(get_ctime()))
    
        



if __name__ == "__main__":
    main()


