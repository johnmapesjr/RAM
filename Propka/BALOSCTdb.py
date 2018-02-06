#!/usr/bin/python3

import glob
import re
import pickle

Xfa = open("../BALOSCTdb.fa").readlines()
P = {}
for i in range(1,len(Xfa),2):
    name = Xfa[i-1].replace('>',"").strip()
    P[name] = {i+1:11.3 for i, x in enumerate(Xfa[i]) if x == 'C'}


for f in glob.glob("Data/*"):
    name = f.replace("Data/","").replace(".B99990001.pka","")
    if name in P:
        f = open(f)
        for line in f:
            if re.match('   CYS', line):
                fields = line.split()
                residue_number = int(fields[1])
                pka = float(fields[3])
                P[name][residue_number] = pka

print(P)

#Calculate Mean
S, n = 0, 0
for k in P:
    for v in P[k].values():
        if v < 99:
            S += v
            n += 1
print(S/n)

pickle.dump(P,open("PropkaBALOSCTdb.dat","wb"))

