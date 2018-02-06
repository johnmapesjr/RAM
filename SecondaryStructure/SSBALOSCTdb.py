#!/usr/bin/python3
import pickle
import glob


AA=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
SS = dict()
n=6
Xfa = open("../BALOSCTdb.fa").readlines()
for i in range(1,len(Xfa),2):
    name = Xfa[i-1].replace('>',"").strip()
    SS[name] = {i+1:[0]*3*(n*2+1) for i, x in enumerate(Xfa[i]) if x == 'C'}

for f in glob.glob("Data/*ss2"):
    name = f.replace("Data/","").replace(".ss2","")
    if name in SS:
        for cys in SS[name]:
            psi = open(f)
            psi.readline()
            psi.readline()
            for line in psi:
                current_residue = int(line.split()[0])
                if cys - n <= current_residue <= cys + n:
                    residue_ss                  = [float(x) for x in line.split()[3:6]]
                    start                       = (n+(cys-current_residue))*3
                    end                         = start+3
                    SS[name][cys][start:end]    = residue_ss

print(SS)
pickle.dump(SS,open("SSBALOSCTdb.dat","wb"))
