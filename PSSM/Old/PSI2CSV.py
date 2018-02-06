#!/usr/bin/python3

from sys import argv, stderr, exit
from math import log, exp
import pickle
import itertools
from hashlib import md5

A={'G':0,'A':1,'V':2,'L':3,'I':4,'P':5,'F':6,'Y':7,'W':8,'S':9,'T':10,'C':11,'M':12,'N':13,'Q':14,'K':15,'R':16,'H':17,'D':18,'E':19}
#working_directory = 'psi_iter3'
working_directory = 'psi_evalue005'
#working_directory = 'psi_trembl_evalue005'

def num_there(s):
    return any(x.isdigit() for x in s)

DIST = pickle.load(open("Modeller100/DIST.dat","rb"))
CSP = pickle.load(open("CSP/CSP.dat","rb"))
#PSS = pickle.load(open("PSIPRED/PSS.dat","rb"))
input_file=open("sp39converted.txt", 'r')
k=6
bonds=[]
bonding=[]
for line in input_file:
    if line == '\n':
        bonds=[]
        bonding=[]
    elif num_there(line):
        line = line.split()
        try:
            a=int(line[0])
            b=int(line[1])
        except:
            a=0
            b=0
        bonds.append((a,b))
        bonding.append(a)
        bonding.append(b)
    else:
        fasta=line.replace("\n","")
        filename = md5(fasta.encode('utf-8')).hexdigest()
        print(filename)
        P = pickle.load( open('%s/%s.dat' % (working_directory, filename), 'rb' ) )
        output_file = open("%s/%s.csv" % (working_directory, filename), 'w')
        if len(P) == 0: continue
        bond_pairs = itertools.combinations(bonding,2)
        bondset = [set(x) for x in bonds]
        for bond in bond_pairs:
            if not all(c in A for c in fasta): continue
            if (bond[0]-1)<0 or (bond[1]-1)<0 or \
            (bond[0]-1) > len(fasta)-1 or (bond[1]-1) > len(fasta)-1: continue
            if fasta[ bond[0]-1 ] == 'C' and fasta[ bond[1]-1 ] == 'C':
                #uncomment to enable pss , did not improve score although was 60% stand alone
#                data=[0]*((k*4+2)*20+(k*4+2)*3+3)
#                data[((k*4+2)*20):((k*4+2)*20+(k*4+2)*3)] = PSS[filename][tuple(sorted((bond[0]-1, bond[1]-1)))] #PSS
#                data[((k*4+2)*20+(k*4+2)*3)] = abs(bond[0]-1 - bond[1]-1) #DOC
#                #Modeller
#                if filename in DIST:
#                    mod_dist = DIST[filename][tuple(sorted((bond[0],bond[1])))]
#                    data[((k*4+2)*20)+(k*4+2)*3+1]= mod_dist
#                else:
#                     data[((k*4+2)*20+(k*4+2)*3)+1] = 0
#                data[((k*4+2)*20+(k*4+2)*3)+2] = CSP[filename][tuple(sorted((bond[0]-1, bond[1]-1)))] #CSP
                data=[0]*((k*4+2)*20+3)
                data[((k*4+2)*20)] = abs(bond[0]-1 - bond[1]-1) #DOC
                #Modeller
                if filename in DIST:
                    mod_dist = DIST[filename][tuple(sorted((bond[0],bond[1])))]
                    data[((k*4+2)*20)+1]= mod_dist
                else:
                     data[((k*4+2)*20+1)] = 0
                data[((k*4+2)*20+2)] = CSP[filename][tuple(sorted((bond[0]-1, bond[1]-1)))] #CSP
                aacount=[1]*(k*4+2)
                #for psi in P:
                for p,psi in enumerate(P,1):
                    left1 = psi[ ((bond[0]-1)-k):(bond[0]-1)]
                    right1 = psi[ ((bond[0]-1)+1):(bond[0]-1)+1+k]
                    left2 = psi[ ((bond[1]-1)-k):(bond[1]-1)]
                    right2 = psi[ ((bond[1]-1)+1):(bond[1]-1)+1+k]
                    #- pad left
                    left1 = ['-']*(k - len(left1)) + left1
                    left2 = ['-']*(k - len(left2)) + left2
                    #- pad right
                    right1 = right1 + ['-']*(k - len(right1))
                    right2 = right2 + ['-']*(k - len(right2))
                    try:
                        cysteine1 = left1 + list(psi[ (bond[0] - 1) ]) + right1
                    except:
                        cysteine1 = (k*2+1)*['-']
                    try:
                        cysteine2 = left2 + list(psi[ (bond[1] - 1) ]) + right2
                    except:
                        cysteine2 = (k*2+1)*['-']
                    if len(cysteine1) == (k*2+1) and len(cysteine2) == (k*2+1):
                        match=cysteine1+cysteine2
                        for i,c in enumerate(match):
                            if c not in A: continue
                            aacount[i] = aacount[i] + 1
                            for j in range(20):
                                if j==A[c]:
                                    data[i*20+j] = data[i*20+j] + 1 * p
                for i in range((k*4+2)*20):
                    #data[i] = data[i] / len(P)
                    data[i] = data[i] / aacount[i//20] 
                    #if (data[i] == 0):
                    #    data[i] = -9999
                    #else:
                    #    data[i] = log(data[i]/0.05,2)
                if set(bond) in bondset:
                    output_file.write("1,"+",".join(str(x) for x in data)+"\n")
                else:
                    output_file.write("0,"+",".join(str(x) for x in data)+"\n")
