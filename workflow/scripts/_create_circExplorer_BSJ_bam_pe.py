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
Logic (for PE reads):
Each BSJ is represented by a 3 alignments in the output BAM file.
Alignment 1 is complete alignment of one of the reads in pair and 
Alignments 2 and 3 are split alignment of the mate at two distinct loci on the same reference 
chromosome.
These alignments are grouped together by the "HI" tags in SAM file. For example, all 3 
alignments for the same BSJ will have the same "HI" value... something like "HI:i:1".
BSJ alignment sam bitflag combinations can have 8 different possibilities, 4 from sense strand
and 4 from anti-sense strand:
1. 83,163,2129
2. 339,419,2385
#           R1.2						      R1.1
#         |<------							<------|
#     5'--|----------------------------------------|-----3'
#     3'--|----------------------------------------|-----5'
#         |                 ------>				   |
#         |                    R2   			   |
#         |                         			   |
#         |<------------------BSJ----------------->|
3. 83,163,2209
4. 339,419,2465
#         						  R1									  
#       						<------									
#     5'--|------------------------------------------|---3'
#     3'--|------------------------------------------|---5'
#         |------>							  ------>|
#         | R2.2								R2.1 | 
#         |                                          |
#         |<-----------------BSJ-------------------->|
5. 99,147,2193
6. 355, 403, 2449
#           R2.1						      R2.2
#         |<------							<------|
#     5'--|----------------------------------------|-----3'
#     3'--|----------------------------------------|-----5'
#         |                 ------>				   |
#         |                    R1   			   |
#         |                         			   |
#         |<------------------BSJ----------------->|
7. 99,147,2145
8. 355, 403, 2401
#         						  R2									  
#       						<------									
#     5'--|------------------------------------------|---3'
#     3'--|------------------------------------------|---5'
#         |------>							  ------>|
#         | R1.2								R1.1 | 
#         |                                          |
#         |<-----------------BSJ-------------------->|
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
        self.bitflags=list()
        self.bitid=""
        self.strand="."
        self.start=-1
        self.end=-1
        self.refcoordinates=dict()
        self.isread1=dict()
        self.isreverse=dict()
        self.issecondary=dict()
        self.issupplementary=dict()
    
    def __str__(self):
        s = "readid: %s"%(self.readid)
        s = "%s\tbitflags: %s"%(s,self.bitflags)
        s = "%s\tbitid: %s"%(s,self.bitid)
        for bf in self.bitflags:
            s = "%s\t%s\trefcoordinates: %s"%(s,bf,", ".join(list(map(lambda x:str(x),self.refcoordinates[bf]))))
        return s

    def set_refcoordinates(self,bitflag,refpos):
        self.refcoordinates[bitflag]=refpos
    
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
    
    def append_alignment(self,read):
        self.alignments.append(read)
    
    def append_bitflag(self,bf):
        self.bitflags.append(bf)
    
    # def extend_ref_positions(self,refcoords):
    # 	self.refcoordinates.extend(refcoords)
    
    def generate_bitid(self):
        bitlist=sorted(self.bitflags)
        self.bitid="##".join(list(map(lambda x:str(x),bitlist)))
# 		self.bitid=str(bitlist[0])+"##"+str(bitlist[1])+"##"+str(bitlist[2])
    
    def get_strand(self):
        if self.bitid=="83##163##2129":
            self.strand="+"
        elif self.bitid=="339##419##2385":
            self.strand="+"
        elif self.bitid=="83##163##2209":
            self.strand="+"
        elif self.bitid=="339##419##2465":
            self.strand="+"		
        elif self.bitid=="99##147##2193":
            self.strand="-"
        elif self.bitid=="355##403##2449":
            self.strand="-"
        elif self.bitid=="99##147##2145":
            self.strand="-"
        elif self.bitid=="355##403##2401":
            self.strand="-"
        elif self.bitid=="16##2064":
            self.strand="+"
        elif self.bitid=="272##2320":
            self.strand="+"
        elif self.bitid=="0##2048":
            self.strand="-"
        elif self.bitid=="256##2304":
            self.strand="-"
        elif self.bitid=="153##2201":
            self.strand="-"
        else:
            self.strand="."
    
    def flip_strand(self):
        if self.strand=="+":self.strand="-"
        if self.strand=="-":self.strand="+"

    def validate_BSJ_read(self,junctions):
        """
        Checks if read is truly a BSJ originitor.
        * Defines left, right and middle alignments
        * Left and right alignments should not overlap
        * Middle alignment should be between left and right alignments
        """
        if len(self.bitid.split("##"))==3:
            left=-1
            right=-1
            middle=-1
            if self.bitid=="83##163##2129":
                left=2129
                right=83
                middle=163
            if self.bitid=="339##419##2385":
                left=2385
                right=339
                middle=419				
            if self.bitid=="83##163##2209":
                left=163
                right=2209
                middle=83
            if self.bitid=="339##419##2465":
                left=419
                right=2465
                middle=339
            if self.bitid=="99##147##2145":
                left=99
                right=2145
                middle=147
            if self.bitid=="355##403##2401":
                left=355
                right=2401
                middle=403
            if self.bitid=="99##147##2193":
                left=2193
                right=147
                middle=99
            if self.bitid=="355##403##2449":
                left=2449
                right=403
                middle=355
            # print(left,right,middle)
            if left == -1 or right == -1 or middle == -1:
                return False
            chrom = self.refname
            # print("validate_BSJ_read",self.readid,left,right,middle)
            # print("validate_BSJ_read",self.readid,self.refcoordinates[left][0],self.refcoordinates[left][-1])
            # print("validate_BSJ_read",self.readid,self.refcoordinates[right][0],self.refcoordinates[right][-1])
            # print("validate_BSJ_read",self.readid,self.refcoordinates[middle][0],self.refcoordinates[middle][-1])
            leftmost = str(self.refcoordinates[left][0])
            rightmost = str(self.refcoordinates[right][-1])
            possiblejid = chrom+"##"+leftmost+"##"+rightmost+"##"+self.strand
            # print("validate_BSJ_read",self.readid,possiblejid)
            if possiblejid in junctions:
                self.start = leftmost
                self.end = str(int(rightmost) + 1)    # this will be added to the BED file
                return True
        else:
            return False
            
    
    
    def get_bsjid(self):
        t=[]
        t.append(self.refname)
        t.append(self.start)
        t.append(self.end)
        t.append(self.strand)
        return "##".join(t)
    
    def write_out_reads(self,outbam):
        for r in self.alignments:
            outbam.write(r)
        
            
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
        description="""Extracts PE BSJs from STAR2p output Chimeric BAM file. It also adds
        unique read group IDs to each read. This RID is of the format <chrom>##<start>##<end>
        where the chrom, start and end represent the BSJ the read is depicting.
        """
    )
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=str,
        help="Input Chimeric-only STAR2p BAM file")
    parser.add_argument('-t','--sample_counts_table', dest='countstable', type=str, required=True,
        help='circExplore per-sample counts table')	# get coordinates of the circRNA
    parser.add_argument("-s",'--sample_name', dest='samplename', type=str, required=False, default = 'sample1',
        help='Sample Name: SM for RG')
    parser.add_argument("-l",'--library', dest='library', type=str, required=False, default = 'lib1',
        help='Sample Name: LB for RG')
    parser.add_argument("-f",'--platform', dest='platform', type=str, required=False, default = 'illumina',
        help='Sample Name: PL for RG')
    parser.add_argument("-u",'--unit', dest='unit', type=str, required=False, default = 'unit1',
        help='Sample Name: PU for RG')
    parser.add_argument("-o","--outbam",dest="outbam",required=True,type=argparse.FileType('w'),
        help="Output bam file ... both strands")
    parser.add_argument("-p","--plusbam",dest="plusbam",required=True,type=argparse.FileType('w'),
        help="Output plus strand bam file")
    parser.add_argument("-m","--minusbam",dest="minusbam",required=True,type=argparse.FileType('w'),
        help="Output plus strand bam file")
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
    # print("Opening...")
    # print(args.inbam)
    print("%s | Reading...alignments!..."%(get_ctime()))
    count=0
    count2=0
    for read in samfile.fetch():
        count+=1
        if debug: print(read,read.reference_id,read.next_reference_id)    
        if read.reference_id != read.next_reference_id: continue    # only works for PE ... for SE read.next_reference_id is -1
        count2+=1
        rid=get_uniq_readid(read)                           # add the HI number to the readid
        if debug:print(rid)
        if not rid in bigdict:
            bigdict[rid]=Readinfo(rid,read.reference_name)
        # bigdict[rid].append_alignment(read)                 # since rid has HI number included ... this separates alignment by HI
        bitflag=get_bitflag(read)
        if debug:print(bitflag)
        bigdict[rid].append_bitflag(bitflag)                # each rid can have upto 3 lines in the BAM with each having its own bitflag ... collect all bigflags in a list here 
        refpos=list(filter(lambda x:x!=None,read.get_reference_positions(full_length=True)))
        bigdict[rid].set_refcoordinates(bitflag,refpos)     # maintain a list of reference coordinated that are "aligned" for each bitflag in each rid alignment
        # bigdict[rid].set_read1_reverse_secondary_supplementary(bitflag,read)
        if debug:print(bigdict[rid])
    print("%s | Done reading %d chimeric alignments. [%d same chrom chimeras]"%(get_ctime(),count,count2))
    samfile.reset()
    print("%s | Writing BAMs"%(get_ctime()))
    print("%s | Re-Reading...alignments!..."%(get_ctime()))
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
        if read.reference_id != read.next_reference_id: continue
        rid=get_uniq_readid(read)
        if rid in bigdict:
            bigdict[rid].generate_bitid()                       # separate all bitflags for the same rid with ## and create a unique single bitflag ... bitflags are pre-sorted
            if debug:print(bigdict[rid])                        
            bigdict[rid].get_strand()                           # use the unique aggregated bitid to extract the strand information ... all possible cases are explicitly covered
            bigdict[rid].flip_strand()                          # strands are flipped than those reported in the counts table .. hence flipping!
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


