import argparse
from os.path import join

def _get_counts_file_path(sampledir,s,prog):
    if prog=="circExplorer":
        return join(sampledir,"circExplorer",s+".circExplorer.counts_table.tsv")
    elif prog=="ciri":
        return join(sampledir,"ciri",s+".ciri.out.filtered")
    elif prog=="dcc":
        return join(sampledir,"DCC",s+".dcc.counts_table.tsv.filtered")
    elif prog=="mapsplice":
        return join(sampledir,"MapSplice",s+".mapsplice.counts_table.tsv.filtered")
    elif prog=="nclscan":
        return join(sampledir,"NCLscan",s+".nclscan.counts_table.tsv.filtered")
    elif prog=="circrnafinder":
        return join(sampledir,"circRNA_finder",s+".circRNA_finder.counts_table.tsv.filtered")


def main() :
    parser = argparse.ArgumentParser(description='Merge per sample Counts from different circRNA detection tools')
    parser.add_argument('--pyscript', dest='pyscript', type=str, required=True,
        help='python script to be run')
    parser.add_argument('--sampledir', dest='sampledir', type=str, required=True,
        help='generally <workdir>/<results>/<sampleName>')
    parser.add_argument('--dcc', dest='dcc', type=int, required=False, default=0,
        help='n_run_dcc')
    parser.add_argument('--mapsplice', dest='mapsplice', type=int, required=False,default=0,
        help='n_run_mapslice')
    parser.add_argument('--nclscan', dest='nclscan', type=int, required=False,default=0,
        help='n_run_nclscan')
    parser.add_argument('--circrnafinder', dest='circrnafinder', type=int, required=False, default=0,
        help='n_run_circrnafinder')
    parser.add_argument('--samplename', dest='samplename', type=str, required=True,
        help='Sample Name')
    parser.add_argument('--min_read_count_reqd', dest='minreads', type=int, required=False, default=2,
        help='Read count threshold..circRNA with lower than this number of read support are excluded! (default=2)')
    parser.add_argument("--reffa",dest="reffa",required=True,type=str,
        help="reference fasta file path")
    parser.add_argument('--pyscriptoutfile',required=True,help='merged table')
    parser.add_argument('--outscript',required=True,help='output bash script')
    args = parser.parse_args()
    
    sd=args.sampledir
    sn=args.samplename
    parameters=""
    parameters+=" --circExplorer "+_get_counts_file_path(sd,sn,'circExplorer')
    parameters+=" --ciri "+_get_counts_file_path(sd,sn,'ciri')
    if args.dcc==1: parameters+=" --dcc "+_get_counts_file_path(sd,sn,'dcc')
    if args.mapsplice==1: parameters+=" --mapsplice "+_get_counts_file_path(sd,sn,'mapsplice')
    if args.nclscan==1: parameters+=" --nclscan "+_get_counts_file_path(sd,sn,'nclscan')
    if args.circrnafinder==1: parameters+=" --circrnafinder "+_get_counts_file_path(sd,sn,'circrnafinder')
    parameters+=" --reffa "+args.reffa
    parameters+=" --min_read_count_reqd "+str(args.minreads)
    parameters+=" --samplename "+sn
    parameters+=" -o "+args.pyscriptoutfile
    
    with open(args.outscript,'w') as outscript:
        outscript.write("python3 "+args.pyscript+parameters)
    outscript.close()

if __name__ == "__main__":
    main()