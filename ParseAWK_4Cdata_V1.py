#!/usr/bin/env python

import sys,re

filename = sys.argv[1]

with open(filename,"r") as f:
    for s in f:
        sr = s.strip().split('\t')
        if 'IsValid' in sr[2] and not any([w in sr[8] for w in ['_and_','BothRepeats','notValid']]):
            patt = re.search(r'([^:]+):(\d+)-(\d+)',sr[1])
            if patt:
                coord = patt.groups()
#                if float(sr[11])>0.0:
                #yield (coord[0], int(coord[1])-1, int(coord[2]), float(sr[11]))
                print coord[0]+"\t"+str(int(coord[1])-1)+"\t"+str(int(coord[2]))+"\t"+str(float(sr[11]))
