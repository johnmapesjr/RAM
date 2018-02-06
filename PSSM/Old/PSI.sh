#!/bin/bash
#working_directory="psi_evalue005"
working_directory="psi_pssm"
#working_directory="psi_iter3"
#working_directory="psi__trembl_evalue005"
grep '[A-Z]' sp39converted.txt | while read line;do 
    md5=`echo -n $line|md5sum|cut -f1 -d' '`
    tmp=`mktemp`
    echo -n $line > $tmp
    if [ ! -f "$working_directory/$md5" ]; then
        psiblast -db ../PSIBLAST/swissprot.fa -query $tmp -evalue=0.005 -num_threads 4 -out_ascii_pssm "$working_directory/$md5" -num_iterations 2
        #psiblast -db ../PSIBLAST/swissprot.fa -query $tmp -evalue=0.005 -num_threads 4 > "psi_evalue005/$md5"
        #psiblast -db ~/Desktop/cys/PSIBlast/trembl.fa -query $tmp -evalue=0.005 -num_threads 4 > "psi_trembl_evalue005/$md5"
        #psiblast -db ~/Desktop/cys/PSIBlast/swissprot.fa -query $tmp -num_iterations 3 -num_threads 4 -evalue 0.005 > "psi_iter3/$md5"
    fi
    rm $tmp
done 
