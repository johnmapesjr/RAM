#!/usr/bin/python3
import pickle
import glob


AA=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
PSSM = dict()
n=6
Xfa = open("../RSC758.fa").readlines()
for i in range(1,len(Xfa),2):
    name = Xfa[i-1].replace('>',"").strip()
    PSSM[name] = {i+1:[0]*20*(n*2+1) for i, x in enumerate(Xfa[i]) if x == 'C'}

for f in glob.glob("Data/*"):
    name = f.replace("Data/","").replace(".pssm","")
    if name in PSSM:
        for cys in PSSM[name]:
            psi = open(f)
            psi.readline()
            psi.readline()
            psi.readline()
            for line in psi:
                if line == '\n':
                    break
                current_residue = int(line.split()[0])
                if cys - n <= current_residue <= cys + n:
                    residue_pssm                = [int(x) for x in line.split()[2:22]]
                    start                       = (n+(cys-current_residue))*20
                    end                         = start+20
                    PSSM[name][cys][start:end]  = residue_pssm

print(PSSM)
pickle.dump(PSSM,open("PSSMRSC758.dat","wb"))
