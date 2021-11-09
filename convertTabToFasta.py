

# usage: python3 convertTabToFasta.py


# 1. convert tab to fasta

from Bio import SeqIO
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

records = SeqIO.parse(infile, "tab")
count = SeqIO.write(records, outfile, "fasta")
print("Converted %i records" % count)
