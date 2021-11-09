

# usage: python3 convertFastaToTab.py


# 1. convert fasta to tab

from Bio import SeqIO
import sys


#infile="./BLASTtest/getorfPB.10054.fasta"
#outfile="./data/getorfPB.10054.tab"
infile = sys.argv[1]
outfile = sys.argv[2]

records = SeqIO.parse(infile, "fasta")
count = SeqIO.write(records, outfile, "tab")
print("Converted %i records" % count)
