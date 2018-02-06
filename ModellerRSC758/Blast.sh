#!/bin/bash

cat "../RSC758.fa"| while read name; read fasta; do
    name=`echo $name|sed -e 's/>//'`
    echo $name
    echo $fasta
    echo
    t=`mktemp`
    echo $fasta > $t
    psiblast -db ../../PDB/PSIBLAST/pdbaa -query $t -out Data/$name -evalue=100000 -num_iterations=100 -num_threads=4
    rm $t
done
