import argparse
import pandas

# pandas.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Test')
parser.add_argument('-o',dest='outfile',required=True,type=argparse.FileType('w'),help='output table')
parser.add_argument('-i',dest='infile', type=argparse.FileType('r'), nargs='+',required=True,help='input files')
args = parser.parse_args()

print(args.infile)
print(args.outfile)

for f in args.infile:
	intable=pandas.read_csv(f,sep="\t",header=0)
	print(intable.head())

