import argparse
import sys
import HTSeq
import gzip
import copy


class BSJ:
    def __init__(self,linestr):
        l=linestr.strip().split("\t")
        self.chrom=l[0]
        self.start=l[1]
        self.end=l[2]
        self.name=l[3]
        self.score=l[4]
        self.strand=l[5]
        self.bitids=l[6]
        self.rids=l[7]
        self.splice_site_flank_5="" #donor
        self.splice_site_flank_3="" #acceptor
    
    def get_jid(self):
        jid=self.chrom+"##"+str(self.start)+"##"+str(self.end)
        return jid
    
    def add_flanks(self,sequences):
        if self.strand == '+':
            coord = int(self.end)
            self.splice_site_flank_5 = sequences[self.chrom][coord:coord+2]
            coord = int(self.start)
            self.splice_site_flank_3 = sequences[self.chrom][coord-2:coord]
        elif self.strand == '-':
            coord = int(self.end)
            myseq = HTSeq.Sequence(bytes(sequences[self.chrom][coord:coord+2],'utf-8'),"myseq")
            revcomp = myseq.get_reverse_complement().seq.decode('utf-8')
            self.splice_site_flank_3 = revcomp
            coord = int(self.start)
            myseq = HTSeq.Sequence(bytes(sequences[self.chrom][coord-2:coord],'utf-8'),"myseq")
            revcomp = myseq.get_reverse_complement().seq.decode('utf-8')
            self.splice_site_flank_5 = revcomp

    def write_out_BSJ(self,outbed):
        t=[]
        t.append(self.chrom)
        t.append(str(self.start))
        t.append(str(self.end))
        t.append(self.name)
        t.append(str(self.score))
        t.append(self.strand)
        t.append(self.bitids)
        t.append(self.rids)
        t.append("##".join([self.splice_site_flank_5,self.splice_site_flank_3]))
        outbed.write("\t".join(t)+"\n")	        

def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
        description="Append the BSJ Donor##Acceptor column to BSJ bed file. Input BSJ bed file is output from _create_circExplorer_BSJ_bam_pe or _create_circExplorer_BSJ_bam_se scripts."
    )
    parser.add_argument("--reffa",dest="reffa",required=True,type=argparse.FileType('r'),default=sys.stdin,
        help="reference fasta file")
    parser.add_argument("--inbsjbedgz",dest="inbsjbedgz",required=True,type=str,
        help="BSJ BED in gzip format")
    parser.add_argument("--outbsjbedgz",dest="outbsjbedgz",required=True,type=str,
        help="BSJ BED in gzip format")
    args = parser.parse_args()

    print("Reading...reference sequences...")
    sequences = dict((s[1], s[0]) for s in HTSeq.FastaReader(args.reffa, raw_iterator=True))
    print("Done reading...%d sequences!"%(len(sequences)))

    print("Reading/Writing...BSJs...")
    bsjs = dict()
    with gzip.open(args.outbsjbedgz,'wt') as bsjfile:
        with gzip.open(args.inbsjbedgz,'rt') as tfile:
            for l in tfile:
                bsj = BSJ(l)
                bsj.add_flanks(sequences)
                bsj.write_out_BSJ(bsjfile)
                bsjs[bsj.get_jid()]=1

    tfile.close()
    bsjfile.close()
    print("Done reading/writing...%d BSJs!"%(len(bsjs)))
    print("Finished!")

if __name__ == "__main__":
    main()