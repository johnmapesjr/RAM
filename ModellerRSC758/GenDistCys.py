#!/usr/bin/python3

from hashlib import md5
import os.path
import sys
import traceback
import itertools
import numpy as np
import pickle
import re

fa      = open("../RSC758.fa").readlines()
fi      = open("../RSC758.txt").readlines()[1:]
names   = [x.split()[0] for x in fi]
idxs    = [int(x.split()[1]) for x in fi]
DIST    = {(name,idx):[] for name,idx in zip(names,idxs)}
for i in range(1,len(fa),2):
    CYS = dict()
    fasta = fa[i].strip()
    filename = fa[i-1].replace(">","").strip()
    cysteines = [i+1 for i,x in enumerate(list(fasta)) if x == "C"]
    print("processing", filename)
    if (not os.path.isfile('Data/%s.B99990001.pdb' % filename)):
        print('PDB DOESNT EXIST %s' % filename)
        continue
    pdb = open("Data/%s.B99990001.pdb" % filename)
    for line in pdb:
        fields = line.split()
        if "ATOM" == fields[0]:
            if "SG" == fields[2] and "CYS" == fields[3]:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                CYS[int(fields[4])] = [x,y,z]
    for cys_check in cysteines:
        if (filename,cys_check) in DIST:
            for cys in cysteines:
                a       = np.array(CYS[cys_check])
                b       = np.array(CYS[cys])
                dist    = np.linalg.norm(a-b)
                DIST[(filename,cys_check)].append( dist )

picklefile = open("DIST.dat","wb")
pickle.dump(DIST,picklefile)
