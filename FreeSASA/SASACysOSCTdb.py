#!/usr/bin/python3

from hashlib import md5
import os.path
import sys
import traceback
import itertools
import numpy as np
import pickle
import re

fa      = open("../BALOSCTdb.fa").readlines()
fi      = open("../BALOSCTdb.txt").readlines()[1:]
names   = [x.split()[0] for x in fi]
idxs    = [int(x.split()[1]) for x in fi]
SASA    = {(name,idx):9.9 for name,idx in zip(names,idxs)}
for i in range(1,len(fa),2):
    CYS = dict()
    fasta = fa[i].strip()
    filename = fa[i-1].replace(">","").strip()
    cysteines = [i+1 for i,x in enumerate(list(fasta)) if x == "C"]
    print("processing", filename)
    if (not os.path.isfile('Data/%s.B99990001.pdb.sasa' % filename)):
        print('PDB DOESNT EXIST %s' % filename)
        continue
    pdb = open("Data/%s.B99990001.pdb.sasa" % filename)
    for line in pdb:
        fields = line.split()
        if "ATOM" == fields[0]:
            if "SG" == fields[2] and "CYS" == fields[3]:
                CYS[int(fields[4])] = float(fields[-1])
    for cys in cysteines:
        if CYS == {}:
            break
        sasa    = CYS[cys]
        SASA[(filename,cys)] = sasa

picklefile = open("SASAOSCTdb.dat","wb")
pickle.dump(SASA,picklefile)
print(SASA)
#mean
S, n = 0, 0
for v in SASA.values():
    S += v
    n += 1
print(S/n)
print(len(SASA))
