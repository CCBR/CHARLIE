import pysam
import argparse
import os
import gzip
import pprint
import time

def get_ctime():
    return time.ctime(time.time())

pp = pprint.PrettyPrinter(indent=4)


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

def _convertjid(jid):
    jid = jid.split("##")
    chrom = jid[0]
    start = jid[1]
    end = jid[2]
    strand = jid[3]
    read_strand = jid[4]
    strand_info = "."
    if strand==read_strand: strand_info="SS"
    if (strand=="+" and read_strand=="-") or (strand=="-" and read_strand=="+"): strand_info="OS"
    return "##".join([chrom,start,end,strand,strand_info])

def _get_shortjid(jid):
    jid = jid.split("##")
    chrom = jid[0]
    start = jid[1]
    end = jid[2]
    strand = jid[3]
    read_strand = jid[4]
    strand_info = "."
    return "##".join([chrom,start,end,strand])

def _get_jinfo(jid):
    jid = jid.split("##")
    chrom = jid[0]
    start = jid[1]
    end = jid[2]
    strand = jid[3]
    read_strand = jid[4]
    strand_info = "."
    if strand==read_strand: strand_info="SS"
    if (strand=="+" and read_strand=="-") or (strand=="-" and read_strand=="+"): strand_info="OS"
    short_jid = "##".join([chrom,start,end,strand])
    converted_jid = "##".join([chrom,start,end,strand,strand_info])
    return chrom,start,end,strand_info,short_jid,converted_jid,read_strand    

class JID:
    def __init__(self,chrom,start,end,strand):
        self.chrom=chrom
        self.start=start
        self.end=end
        self.strand=strand
        self.ss_linear_count=0
        self.os_linear_count=0
        self.ss_linear_spliced_count=0
        self.os_linear_spliced_count=0

    def increment_linear(self,strand_info):
        if strand_info=="SS": self.ss_linear_count+=1
        if strand_info=="OS": self.os_linear_count+=1

    def increment_linear_spliced(self,strand_info):
        if strand_info=="SS": self.ss_linear_spliced_count+=1
        if strand_info=="OS": self.os_linear_spliced_count+=1


def main():
    # debug = True
    debug = False
    parser = argparse.ArgumentParser(
    )
    # INPUTs
    parser.add_argument("-i","--inbam",dest="inbam",required=True,type=str,
        help="Input BAM file")
    parser.add_argument('-r',"--rid2jid",dest="rid2jid",required=True,type=str,
        help="readID to junctionID lookup")
    parser.add_argument('-t','--sample_counts_table', dest='countstable', type=str, required=True,
        help='circExplore per-sample counts table')	# get coordinates of the circRNA
    parser.add_argument("-s",'--sample_name', dest='samplename', type=str, required=False, default = 'sample1',
        help='Sample Name: SM for RG')
    parser.add_argument('-p',"--pe",dest="pe",required=False,action='store_true', default=False,
        help="set this if BAM is paired end")
    parser.add_argument("-l",'--library', dest='library', type=str, required=False, default = 'lib1',
        help='Sample Name: LB for RG')
    parser.add_argument("-f",'--platform', dest='platform', type=str, required=False, default = 'illumina',
        help='Sample Name: PL for RG')
    parser.add_argument("-u",'--unit', dest='unit', type=str, required=False, default = 'unit1',
        help='Sample Name: PU for RG')
    parser.add_argument('--regions', dest='regions', type=str, required=True,
        help='regions file eg. ref.fa.regions')
    parser.add_argument('--host', dest='host', type=str, required=True,
        help='host name eg.hg38... single value')
    parser.add_argument('--additives', dest='additives', type=str, required=True,
        help='additive name(s) eg.ERCC... comma-separated list... all BSJs in this region are filtered out')
    parser.add_argument('--viruses', dest='viruses', type=str, required=True,
        help='virus name(s) eg.NC_009333.1... comma-separated list')
    # OUTPUTs
    parser.add_argument("-o","--outbam",dest="outbam",required=True,type=str,
        help="Output \"primary alignment near BSJ\" only BAM file")
    parser.add_argument("--outplusbam",dest="outplusbam",required=True,type=str,
        help="Output \"primary alignment near BSJ\" only plus strand BAM file")
    parser.add_argument("--outminusbam",dest="outminusbam",required=True,type=str,
        help="Output \"primary alignment near BSJ\" only minus strand BAM file")
    parser.add_argument("--splicedbam",dest="splicedbam",required=True,type=str,
        help="Output \"primary spliced alignment\" only BAM file")
    parser.add_argument("--splicedbsjbam",dest="splicedbsjbam",required=True,type=str,
        help="Output \"primary spliced alignment near BSJ\" only BAM file")
    parser.add_argument("--splicedbsjplusbam",dest="splicedbsjplusbam",required=True,type=str,
        help="Output \"primary spliced alignment near BSJ\" only plus strand BAM file")
    parser.add_argument("--splicedbsjminusbam",dest="splicedbsjminusbam",required=True,type=str,
        help="Output \"primary spliced alignment near BSJ\" only minus strand BAM file")
    parser.add_argument("--outputhostbams",dest="outputhostbams",required=False,action='store_true', default=False,
        help="Output individual host BAM files")
    parser.add_argument("--outputvirusbams",dest="outputvirusbams",required=False,action='store_true', default=False,
        help="Output individual virus BAM files")
    parser.add_argument("--outdir",dest="outdir",required=False,type=str,
        help="Output folder for the individual BAM files (required only if --outputhostbams or --outputvirusbams is used).")
    parser.add_argument("-c","--countsfound",dest="countsfound",required=True,type=argparse.FileType('w', encoding='UTF-8'),
        help="Output TSV file with counts of junctions found")


    args = parser.parse_args()
    print("%s | Reading...rid2jid!..."%(get_ctime()))
    rid2jid = dict()
    with gzip.open(args.rid2jid,'rt') as tfile:
        for l in tfile:
            l=l.strip().split("\t")
            rid2jid[l[0]]=l[1]
    tfile.close()
    print("%s | Done reading...%d rid2jid's!"%(get_ctime(),len(rid2jid)))

    samfile = pysam.AlignmentFile(args.inbam, "rb")
    samheader = samfile.header.to_dict()
    samheader['RG']=list()
    junctionsfile = open(args.countstable,'r')
    print("%s | Reading...junctions!..."%(get_ctime()))
    count=0
    junction_counts=dict()
    # splicedbsjjid=dict()
    for l in junctionsfile.readlines():
        count+=1
        if "read_count" in l: continue
        l = l.strip().split("\t")
        chrom = l[0]
        start = l[1]
        end = str(int(l[2])-1)
        strand = l[3]
        short_jid  = chrom+"##"+start+"##"+end+"##"+strand        # create a unique junction ID for each line in the BSJ junction file and make it the dict key ... easy for searching!
        jid1 = short_jid+"##SS"                                   # SS=sample strand ... called BSJ and read are on the same strand
        jid2 = short_jid+"##OS"                                   # OS=opposite strand ... called BSJ and read are on opposite strands
        samheader['RG'].append({'ID':jid1 ,  'LB':args.library, 'PL':args.platform, 'PU':args.unit,'SM':args.samplename})
        samheader['RG'].append({'ID':jid2 ,  'LB':args.library, 'PL':args.platform, 'PU':args.unit,'SM':args.samplename})
        # print(short_jid)
        junction_counts[short_jid] = JID(chrom,start,end,strand)
        # splicedbsjjid[jid] = dict()
    junctionsfile.close()
    # exit()
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
    print("%s | Done reading %d junctions."%(get_ctime(),count))
    
    outbam = pysam.AlignmentFile(args.outbam, "wb", header=samheader)
    outplusbam = pysam.AlignmentFile(args.outplusbam, "wb", header=samheader)
    outminusbam = pysam.AlignmentFile(args.outminusbam, "wb", header=samheader)
    splicedbam = pysam.AlignmentFile(args.splicedbam, "wb", header=samheader)
    splicedbsjbam = pysam.AlignmentFile(args.splicedbsjbam, "wb", header=samheader)
    splicedbsjplusbam = pysam.AlignmentFile(args.splicedbsjplusbam, "wb", header=samheader)
    splicedbsjminusbam = pysam.AlignmentFile(args.splicedbsjminusbam, "wb", header=samheader)
    outputbams = dict()
    if args.outputhostbams:
        for h in hosts:
            outbamname = os.path.join(args.outdir,args.samplename+"."+h+".BSJ.bam")
            outputbams[h] = pysam.AlignmentFile(outbamname, "wb", header = samheader)
    if args.outputvirusbams:
        for v in viruses:
            outbamname = os.path.join(args.outdir,args.samplename+"."+v+".BSJ.bam")
            outputbams[v] = pysam.AlignmentFile(outbamname, "wb", header = samheader)            
    lenoutputbams = len(outputbams)
    # pp.pprint(rid2jid)
    print("%s | Opened output BAMs for writing..."%(get_ctime()))
    spliced=dict() # 1=spliced
    splicedbsj=dict()
    count1=0    # total reads
    count2=0    # total reads near BSJ
    count3=0    # total spliced reads
    count4=0    # total spliced reads near BSJ
    print("Reading alignments...")
    mate_already_counted1=dict()
    mate_already_counted2=dict()
    # mate_already_counted3=dict() # not needed as similar to the "spliced" dict
    # mate_already_counted4=dict() # not needed as similar to "spliced" dict have value 2
    last_printed=-1
    for read in samfile.fetch():
        if args.pe and ( read.reference_id != read.next_reference_id ): continue    # only works for PE ... for SE read.next_reference_id is -1
        if args.pe and ( not read.is_proper_pair ): continue
        if read.is_secondary or read.is_supplementary or read.is_unmapped : continue
        rid=read.query_name
# count read if it has not been counted yet
        if not rid in mate_already_counted1:
            mate_already_counted1[rid]=1
            count1+=1
# find cigar tuple, cigar string and generate a cigar string order
# if 0 is followed by 3 in the "cigarstringorder" value (can happen more than once in multi-spliced reads)
# then the read is spliced
        cigar=read.cigarstring
        cigart=read.cigartuples
        cigart=cigart[list(map(lambda z:z[0],cigart)).index(0):]
        cigarstringorder=""
        for j in range(len(cigart)):
            cigarstringorder+=str(cigart[j][0])
# cigarstringorder can be like 034 or 03034 or 03 or 0303
# check if the rid is already found to be spliced ... if not then check if it is
        if not rid in spliced:
            if "03" in cigarstringorder: # aka read is spliced
                count3+=1
                spliced[rid]=1
# check if the rid exists in the rid2jid lookup table
        if rid in rid2jid:  # does this rid have a corresponding BSJ??
# if rid is in rid2jid lookuptable and it is not previously counted then count it as "linear" read for that BSJ
            if not rid in mate_already_counted2:
                mate_already_counted2[rid]=1
                count2+=1
            jid = rid2jid[rid]
            chrom, jstart, jend, strand_info, short_jid, converted_jid, read_strand = _get_jinfo(jid)
            junction_counts[short_jid].increment_linear(strand_info)
            read.set_tag("RG", converted_jid , value_type="Z")
            outbam.write(read)
            if read_strand=="+": outplusbam.write(read)
            if read_strand=="-": outminusbam.write(read)
            if lenoutputbams != 0:
                regionname=_get_regionname_from_seqname(regions,chrom)
                if regionname in hosts and args.outputhostbams:
                    outputbams[regionname].write(read)
                if regionname in viruses and args.outputvirusbams:
                    outputbams[regionname].write(read)
# check if this rid's .. this alignment is spliced!
# rid could be in spliced but this may be an unspliced mate
            if rid in spliced and "03" in cigarstringorder:
                if not rid in splicedbsj:
# CIGAR has match ... followed by skip ... aka spliced read
# find number of splices
# nsplices is the number of times "03" is found in cigarstringorder
# if nsplices is gt than 1 then we have to get the coordinates of all the matches and 
# try to compare each one with the BSJ coordinates
                    nsplices = cigarstringorder.count("03")
                    if nsplices == 1:
                        start=int(read.reference_start)+int(cigart[0][1])+1
                        end=int(start)+int(cigart[1][1])-1
                        # print(start,end,jstart,jend)
                        if abs(int(start)-int(jstart))<3 or abs(int(end)-int(jend))<3: # include 2,1,0,-1,-2
                            junction_counts[short_jid].increment_linear_spliced(strand_info)
                            splicedbsj[rid]=1  # aka read is spliced and is spliced at BSJ
                            count4+=1
                            # splicedbsjjid[jid][rid]=1
                    else:   # read has multiple splicing events
                        for j in range(len(cigart)-1):
                            if cigart[j][0]==0 and cigart[j+1][0]==3:
                                add_coords = 0
                                for k in range(j+1):
                                    add_coords+=int(cigart[k][1])
                                start=int(read.reference_start)+add_coords+1
                                end=int(start)+int(cigart[j+1][1])-1
                                if abs(int(start)-int(jstart))<3 or abs(int(end)-int(jend))<3: # include 2,1,0,-1,-2
                                    junction_counts[short_jid].increment_linear_spliced(strand_info)
                                    splicedbsj[rid]=1  # aka read is spliced and is spliced at BSJ
                                    count4+=1
                                    # splicedbsjjid[jid][rid]=1
                                    break
        if (count1%100000==0) and (last_printed!=count1):
            last_printed=count1
            print("%s | ...Processed %d reads/readpairs (%d  were spliced! %d linear around BSJ! %d spliced at BSJ)"%(get_ctime(),count1,len(spliced),count2,len(splicedbsj)))
    print("%s | Done processing alignments: %d reads/readpairs (%d  were spliced! %d linear around BSJ! %d spliced at BSJ)"%(get_ctime(),count1,len(spliced),count2,len(splicedbsj)))
    if lenoutputbams != 0:
        for k,v in outputbams.items():
            v.close()
    samfile.reset()
    print("%s | Writing spliced BAMs ..."%(get_ctime()))

    for read in samfile.fetch():
        rid = read.query_name
        if rid in spliced : splicedbam.write(read)
        if rid in splicedbsj : 
            jid = rid2jid[rid]
            # converted_jid = _convertjid(jid)
            chrom, jstart, jend, strand_info, short_jid, converted_jid, read_strand = _get_jinfo(jid)
            read.set_tag("RG", converted_jid ,  value_type="Z") 
            splicedbsjbam.write(read)
            if read_strand=="+": splicedbsjplusbam.write(read)
            if read_strand=="-": splicedbsjminusbam.write(read)

    samfile.close()
    outbam.close()
    outplusbam.close()
    outminusbam.close()
    splicedbam.close()
    splicedbsjbam.close()
    splicedbsjplusbam.close()
    splicedbsjminusbam.close()
    print("%s | Closing all BAMs"%(get_ctime()))
    args.countsfound.write("#chrom\tstart\tend\tstrand\tlinear_BSJ_reads_same_strand\tlinear_spliced_BSJ_reads_same_strand\tlinear_BSJ_reads_opposite_strand\tlinear_spliced_BSJ_reads_opposite_strand\n")
    for short_jid in junction_counts.keys():
        chrom=junction_counts[short_jid].chrom
        start=junction_counts[short_jid].start
        end=int(junction_counts[short_jid].end)+1
        strand=junction_counts[short_jid].strand
        ss_linear_count=junction_counts[short_jid].ss_linear_count
        ss_linear_spliced_count=junction_counts[short_jid].ss_linear_spliced_count
        os_linear_count=junction_counts[short_jid].os_linear_count
        os_linear_spliced_count=junction_counts[short_jid].os_linear_spliced_count
        args.countsfound.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n"%(chrom,str(start),str(end),strand,ss_linear_count,ss_linear_spliced_count,os_linear_count,os_linear_spliced_count))
    args.countsfound.close()
    print("%s | DONE!!"%(get_ctime()))


if __name__ == "__main__":
    main()