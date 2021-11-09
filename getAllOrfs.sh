#!/bin/bash

# GIT_REPODIR="/Users/clarapereira/tools/TransDecoder-TransDecoder-v5.5.0"
# refdir="/Users/clarapereira/Dropbox/Barreto_lab/reference"
#
# #gtfPATH="$HOME/Dropbox/Barreto_lab/final_effort_HSCs_AI/reference"
# gtf=gencode.vM1.annotation.gtf.gz
#
# ${GIT_REPODIR}/util/gtf_to_alignment_gff3.pl ${refdir}/${gtf%.*.*}.KNOWNtranscript.gtf > ${refdir}/${gtf%.*.*}.KNOWNtranscript.gff3
#
# ${GIT_REPODIR}/TransDecoder.LongOrfs -t ${refdir}/gencode.vM27.transcripts.fasta
#
#
# gunzip -c ${refdir}/${gtf} | \
# 	grep "chr" |
# 	grep 'gene_status "KNOWN"' |
# 	awk 'BEGIN{OFS="\t";} $3=="transcript" {print}' > ${refdir}/${gtf%.*.*}.KNOWNtranscript.gtf
#
#
#
# ${GIT_REPODIR}/util/gtf_genome_to_cdna_fasta.pl ${refdir}/${gtf%.*.*}.KNOWNtranscript.gtf ${refdir}/GRCm39.genome.fa > ${refdir}/transcripts.fasta

# Strategy:
# 1. Obtain all the orfs (if using the ORF, then use blastx) OR obtain all the peptides (to use blastp)
# 2. Obtain all the known protein sequences from a database
# 3. Do the blast:
# 3.1. Create tab file from each of the fasta files
# 3.2. Simplify the name of each sequence
# 3.3. Revert back to .fasta
# 3.4. Separate query .fasta file into N .fasta files (as many as the sequences)
# 3.5. Run the blastp /blastx and save the output with identical basename as each of the query fasta files
# 3.6. from the output, retain only the best match
# 3.7. import all the files, merging everything --> the output should have the ORF, the peptide and all the columns of the best match
# 3.8. group by gene and retain only the line with the best match; if there are ties, then choose the longest ORF



dir=/Users/clarapereira/Dropbox/Barreto_lab/ORFs
fastadir=orfs_out/orfipy_gencode.vM27.transcripts.fasta_out

# 1a. obtain all the orfs
# get all ORFs with orfipy (from mouse trasncripts):
orfipy ${refdir}/gencode.vM27.transcripts.fasta --dna orfs.fa --min 10 --max 10000 --procs 4 --table 1 --outdir orfs_out

# 1b. obtain all the peptides
# Translate all the ORFs with orfipy:
orfipy ${refdir}/gencode.vM27.transcripts.fasta --pep orfs_peptides.fa --min 50 --procs 4


# Protein sequences were obtained from uniprop
# the file:
gunzip -c /Users/clarapereira/Dropbox\ \(Personal\)/Barreto_lab/ORFs/reference/UP000000589_10090.fasta.gz | head -100


# 3.1. Create tab file from each of the fasta files
# - orf peptides
python3 convertFastaToTab.py ${dir}/${fastadir}/orfs_peptides.fa ${dir}/${fastadir}/orfs_peptides.tab
# - orf nucleotides
python3 convertFastaToTab.py ${dir}/orfs_out/orfs.fa ${dir}/orfs_out/orfs.tab


# 3.2. Simplify the name of each sequence
# 3.3. Revert back to .fasta
python3 convertTabToFasta.py ${dir}/${fastadir}/orfs_peptides.tab ${dir}/${fastadir}/orfs_peptides.simpleIDs.fasta


# 3.4. Separate query .fasta file into N .fasta files (as many as the sequences)
sh fastaToManyFasta.sh ${dir}/${fastadir}/orfs_peptides.simpleIDs.fasta ${dir}/${fastadir}/fasta
sh fastaToManyFasta.sh ${dir}/${fastadir}/orfs_peptidesTEST.simpleIDs.fasta ${dir}/${fastadir}/fasta_test

# 3.5. Run the blastp /blastx and save the output with identical basename as each of the query fasta files
sh blastpOverFiles.sh \
  "${dir}/orfs_out/orfipy_gencode.vM27.transcripts.fasta_out/fasta_test/" \
  "/Users/clarapereira/Dropbox/Barreto_lab/ORFs/reference/UP000000589_10090.fasta"
#
/usr/local/ncbi/blast/bin/blastp -query ${dir}/orfs_out/orfipy_gencode.vM27.transcripts.fasta_out/orfs_2610203C22Rik-202_random_example.fasta -subject ${dir}/reference/UP000000589_10090.fasta -outfmt 10 -out ${dir}/2610203C22Rik-202.fasta.csv
