cat ../RSC758.fa|while read name; read fasta; do
    tmp=`mktemp`
    pssm=`echo Data/$name.pssm|sed -e 's/>//'`
    echo $fasta > $tmp
    echo "Working on $output"
    psiblast -db Database/swissprot.fa -query $tmp -evalue=0.005 -num_threads 4 -out_ascii_pssm $pssm -num_iterations 2
    rm $tmp
done
