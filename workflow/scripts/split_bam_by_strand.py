import pysam
import sys
import argparse

"""
This script splits a BSJ bam file by strand into:
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
		
	def write_out_BSJ(self,outbed):
		t=[]
		t.append(self.chrom)
		t.append(str(self.start))
		t.append(str(self.end))
		t.append(self.name)
		t.append(str(self.score))
		t.append(self.strand)
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
		self.refcoordinates=list()
	
	def append_alignment(self,read):
		self.alignments.append(read)
	
	def append_bitflag(self,bf):
		self.bitflags.append(bf)
	
	def extend_ref_positions(self,refcoords):
		self.refcoordinates.extend(refcoords)
	
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
		else:
			self.strand="U"
	
	def get_start_end(self):
		refcoords=sorted(self.refcoordinates)
		self.start=str(refcoords[0])
		self.end=str(int(refcoords[-1])+1)
	
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
	parser = argparse.ArgumentParser()
	parser.add_argument("-i","--inbam",dest="inbam",required=True,type=argparse.FileType('r'),
		help="Input bam file")
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
	bigdict=dict()
	for read in samfile.fetch():
		rid=get_uniq_readid(read)
		if not rid in bigdict:
			bigdict[rid]=Readinfo(rid,read.reference_name)
		bigdict[rid].append_alignment(read)
		bigdict[rid].append_bitflag(get_bitflag(read))
		bigdict[rid].extend_ref_positions(read.get_reference_positions(full_length=False))
		
	bsjdict=dict()
	for rid in bigdict.keys():
		bigdict[rid].generate_bitid()
		bigdict[rid].get_strand()
		bigdict[rid].get_start_end()
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
	
	for bsjid in bsjdict.keys():
		bsjdict[bsjid].write_out_BSJ(args.bed)
		
	plusfile.close()
	minusfile.close()
	samfile.close()
	args.bed.close()
	
		



if __name__ == "__main__":
    main()


