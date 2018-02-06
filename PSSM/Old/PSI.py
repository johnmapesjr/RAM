#!/usr/bin/python3
import sys
import pickle
from hashlib import md5

f=open('fastas')
#working_directory = 'psi_iter3'
working_directory = 'psi_evalue005'
#working_directory = 'psi_trembl_evalue005'

for fasta in f:
    fasta = fasta.replace('\n','')
    filename = md5(fasta.encode('utf-8')).hexdigest()
    print(filename)
    psi = open('%s/%s' % (working_directory,filename))
    psilines = psi.readlines()
    P = list()
    P.append(list(fasta))
    j = 0
    for i in range(len(psilines)):
        if "Score =" in psilines[i]:
            P.append(['-']*len(fasta))
            j = j + 1
            k = i + 1
            begin = False
            q=list()
            s=list()
            while( "Score = " not in psilines[k] and "BLOSUM62" not in psilines[k]):
                k = k + 1
                if "Query  " in psilines[k]:
                    if not begin:
                        begin   = int(psilines[k].split()[1])
                    q       += psilines[k].split()[2]
                    end     = int(psilines[k].split()[3])
                    s       += psilines[k+2].split()[2]
            remove  = [x for x, ltr in enumerate(q) if ltr == '-']
            s       = [ltr for x, ltr in enumerate(s) if x not in remove]
            P[j][begin-1:end] = s[0:(end-begin+1)]
            i = k
    psipickle = open(('%s/%s.dat' % (working_directory, filename)),'wb')
    pickle.dump(P,psipickle)
