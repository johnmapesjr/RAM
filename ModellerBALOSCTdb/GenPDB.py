#!/usr/bin/python3

from hashlib import md5
from modeller import *
from modeller.automodel import *    # Load the automodel class
import os.path
import sys
import traceback
import re
#from shutil import copyfile

def mainloop(filename,fasta,cutoff):
    fasta       = fasta.strip()
#    filename    = md5(fasta.encode('utf-8')).hexdigest()
    if (os.path.isfile('Data/%s.B99990001.pdb' % filename)):
        return 'EXISTS %s' % filename
    seg_file    = open("Data/%s.seg" % filename, "w")
    blast       = open("Data/"+filename)
    fasta_len   = len(fasta)
    matches     = set()
    log_file    = open("logpdb.txt","a")
    seg_file.write(">P1;%s\n" % filename)
    seg_file.write("sequence:::::::::\n")
    seg_file.write("%s*\n" % fasta)
    i = 0 
    for line in blast:
        if re.match('^>pdb', line):
            match = line.split('|')[1].split('|')[0]
            match = match.lower()
        if "Identities = " in line:
            match_len   = int(line.split('/')[0].split()[2])
            identity    = match_len / fasta_len
            identity = int(line.split('%')[0].split('(')[1])/100
            if (identity < 1.1 and i < cutoff):
                if (identity == 1 and fasta_len == match_len):
                    log_file.write("----------------------------------------------------------------------------------------\n")
                    log_file.write("----------------------------------------------------------------------------------------\n")
                    log_file.write("SAME PDB\n%s\n%s\n%s\n" % (fasta, match, filename))
                    log_file.write("----------------------------------------------------------------------------------------\n")
                    log_file.write("----------------------------------------------------------------------------------------\n")
                    #include pdb because is available for prediction
#                    copyfile('/root/Desktop/3TB/pdb/%s/pdb%s.ent.gz' % (match[1:3],match),'Data/%s.B99990001.pdb' % filename)
#                    return
                i += 1
                check = '/root/Desktop/3TB/pdb/%s/pdb%s.ent.gz' % (match[1:3],match)
                if match not in matches and os.path.isfile(check):
                    seg_file.write('>P1;%s\n' % match)
                    seg_file.write('structure:%s:FIRST:@:END:::::\n' % match)
                    seg_file.write('*\n')
                    matches.add(match)
    seg_file.close()
    log.verbose()
    env = environ()
    env.io.atom_files_directory = ['/root/Desktop/3TB/pdb/']
    matches = list(matches)
    if (len(matches)==1):
        matches = matches[0]
    a = automodel(env,
                  # file with template codes and target sequence
                  alnfile  = 'Data/%s.seg' % filename,
                  # PDB codes of the templates
                  knowns   = matches,
                  # code of the target
                  sequence = filename)
    a.max_molpdf = 100000000
    a.auto_align()                      # get an automatic alignment
    a.make()                            # do comparative modeling

f=open("../BALOSCTdb.fa").readlines()
for i in range(1,len(f),2):
    fasta = f[i]
    filename = f[i-1].replace(">","").strip() 
    for cutoff in range(2,0,-1):
        try:
            print(mainloop(filename, fasta, cutoff))
            break
        except Exception as e:
            print('Too high cutoff @ %s or %s' % (cutoff, e,))
            traceback.print_exc()


