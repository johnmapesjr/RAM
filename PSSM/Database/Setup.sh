#/bin/bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta swissprot.fa
makeblastdb -in swissprot.fa  -dbtype prot
