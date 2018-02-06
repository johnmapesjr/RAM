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
QUAL    = dict()
for i in range(1,len(fa),2):
    fasta = fa[i].strip()
    filename = fa[i-1].replace(">","").strip()
    with open("Data/%s" % filename, "r") as aln:
        for line in aln:
            if "Length=" in line:
                length = int(line.split("=")[1])
                break
        for line in aln:
            if "Identities = " in line:
                identity = int(line.split("Identities = ")[1].split("/")[0])
                break
    QUAL[filename] = identity/length

picklefile = open("Quality.dat","wb")
pickle.dump(QUAL,picklefile)
#print(QUAL)
for k,v in QUAL.items():
    print(v)
