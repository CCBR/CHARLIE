import pysam
import sys
import argparse

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
        self.bitids=list()
        self.rids=list()
        
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
        self.bitids.append(bitid)

    def append_rid(self,rid):
        self.rids.append(rid)
        
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
        
class Readinfo:
    def __init__(self,readid,rname):
        self.readid=readid
        self.refname=rname
        self.alignments=list()
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
            self.strand="U"

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
            leftmost = str(self.refcoordinates[left][0])
            rightmost = str(self.refcoordinates[right][-1])
            possiblejid = chrom+"##"+leftmost+"##"+rightmost
            if possiblejid in junctions:
                self.start = leftmost
                self.end = str(int(rightmost) + 1)    # this will be added to the BED file
                return True
        else:
            return False
            
    
    # def get_start_end(self):
    #     refcoordinates=self.refcoordinates
    #     isread1=self.isread1
    #     if len(self.isread1)!=3:			
    #         refcoords=[]
    #         for i in refcoordinates.keys():
    #             refcoords.extend(refcoordinates[i])
    #     else:			
    #         l=[]
    #         for i in isread1.keys():
    #             l.append(isread1[i])
    #         Ycount=l.count("Y")
    #         Ncount=l.count("N")
    #         if Ycount>Ncount:
    #             useread1="Y"
    #         else:
    #             useread1="N"
    #         refcoords=[]
    #         for i in refcoordinates.keys():
    #             if isread1[i]==useread1:
    #                 refcoords.extend(refcoordinates[i])
    #     refcoords=sorted(refcoords)
    #     self.start=str(refcoords[0])
    #     self.end=str(int(refcoords[-1])+1)
    
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



def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=argparse.FileType('r'),
        help="Input bam file")
    parser.add_argument('--sample_counts_table', dest='countstable', type=str, required=True,
                    help='circExplore per-sample counts table')	# get coordinates of the circRNA
    parser.add_argument("-p","--plusbam",dest="plusbam",required=True,type=argparse.FileType('w'),
        help="Output plus strand bam file")
    parser.add_argument("-m","--minusbam",dest="minusbam",required=True,type=argparse.FileType('w'),
        help="Output plus strand bam file")
    parser.add_argument("-b","--bed",dest="bed",required=True,type=argparse.FileType('w', encoding='UTF-8'),
        help="Output BSJ bed file (with strand info)")
    args = parser.parse_args()		
    samfile = pysam.AlignmentFile(args.inbam, "rb")
    plusfile = pysam.AlignmentFile(args.plusbam, "wb", template=samfile)
    minusfile = pysam.AlignmentFile(args.minusbam, "wb", template=samfile)
# 	bsjfile = open(args.bed,"w")
    junctionsfile = open(args.countstable,'r')
    junctions=dict()
    print("Reading...junctions!...")
    for l in junctionsfile.readlines():
        if "read_count" in l: continue
        l = l.strip().split("\t")
        chrom = l[0]
        start = l[1]
        end = str(int(l[2])-1)
        jid = chrom+"##"+start+"##"+end                     # create a unique junction ID for each line in the BSJ junction file and make it the dict key ... easy for searching!
        junctions[jid]=1
    junctionsfile.close()
    print("Done reading %d junctions."%(len(junctions)))

    bigdict=dict()
    print("Reading...alignments!...")
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
        bigdict[rid].append_alignment(read)                 # since rid has HI number included ... this separates alignment by HI
        bitflag=get_bitflag(read)
        if debug:print(bitflag)
        bigdict[rid].append_bitflag(bitflag)                # each rid can have upto 3 lines in the BAM with each having its own bitflag ... collect all bigflags in a list here 
        refpos=list(filter(lambda x:x!=None,read.get_reference_positions(full_length=True)))
        bigdict[rid].set_refcoordinates(bitflag,refpos)     # maintain a list of reference coordinated that are "aligned" for each bitflag in each rid alignment
        # bigdict[rid].set_read1_reverse_secondary_supplementary(bitflag,read)
        if debug:print(bigdict[rid])
    print("Done reading %d chimeric alignments. [%d same chrom chimeras]"%(count,count2))

    print("Writing BAMs")
    bsjdict=dict()
    bitid_counts=dict()
    for rid in bigdict.keys():
        bigdict[rid].generate_bitid()                       # separate all bitflags for the same rid with ## and create a unique single bitflag ... bitflags are pre-sorted
        if debug:print(bigdict[rid])                        
        bigdict[rid].get_strand()                           # use the unique aggregated bitid to extract the strand information ... all possible cases are explicitly covered
        if not bigdict[rid].validate_BSJ_read(junctions=junctions): # ensure that the read alignments leftmost and rightmost coordinates match with one of the BSJ junctions... if yes then that rid represents a BSJ. Also add start and end to the BSJ object
            continue
        # bigdict[rid].get_start_end()
        # print(bigdict[rid])
        if bigdict[rid].strand=="+":
            bigdict[rid].write_out_reads(plusfile)
        if bigdict[rid].strand=="-":
            bigdict[rid].write_out_reads(minusfile)
        bsjid=bigdict[rid].get_bsjid()
        if not bsjid in bsjdict:
            bsjdict[bsjid]=BSJ()
            bsjdict[bsjid].set_chrom(bigdict[rid].refname)
            bsjdict[bsjid].set_start(bigdict[rid].start)
            bsjdict[bsjid].set_end(bigdict[rid].end)
            bsjdict[bsjid].set_strand(bigdict[rid].strand)
        bsjdict[bsjid].plusone()
        bsjdict[bsjid].append_bitid(bigdict[rid].bitid)
        if not bigdict[rid].bitid in bitid_counts:
            bitid_counts[bigdict[rid].bitid]=0
        bitid_counts[bigdict[rid].bitid]+=1
        bsjdict[bsjid].append_rid(rid)
    print("Done!")	
    for b in bitid_counts.keys():
        print(b,bitid_counts[b])
    print("Writing BED")
    for bsjid in bsjdict.keys():
        bsjdict[bsjid].write_out_BSJ(args.bed)
        
    plusfile.close()
    minusfile.close()
    samfile.close()
    args.bed.close()
    print("ALL Done!")
    
        



if __name__ == "__main__":
    main()


