#!/bin/bash

infile=${1}
outdir=${2}

mkdir -p ${outdir}

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${outdir}/${line#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < ${infile}
