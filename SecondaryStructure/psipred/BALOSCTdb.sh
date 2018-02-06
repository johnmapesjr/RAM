#!/bin/bash

cat ../../BALOSCTdb.fa| while read name; read fasta; do
    name=`echo $name|sed -e 's/>//'`
    echo $fasta > ../Data/$name.fa
    ./runpsipred ../Data/$name.fa
    rm $name.horiz $name.ss ../Data/*fa
    mv $name.ss2 ../Data/
done
