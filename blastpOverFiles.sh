#!/bin/bash

infile=${1}
db=${2}


for file in "${infile}"/*.fa; do
	[ -e "${file}" ] || continue
	echo ${file}
	echo "..."
	/usr/local/ncbi/blast/bin/blastp \
	  -query ${file} \
	  -subject ${db} \
	  -outfmt 10 \
	  -out ${file}_blastpResult.fasta.tsv
done

echo "Done."
