#!/usr/bin/env python


"""
=========================
To parse the 4C-rmsk library
=========================
scripts adapted from bbcflib.createlib

"""


import os,sys


resfile = sys.argv[2]
outf = open(resfile,'w')


with open(sys.argv[1]) as f:
	for line in f:
		s = line.split('\t')
		s_split = s[3].split('|')
		infos = '|'.join(s_split[0:(len(s_split)-4)]+list(s[4:8]))
		outf.write(s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + infos)




#for s in sorted_stream(coverout.read(),[chrom]):
#	s_split = s[3].split('|')
#	infos = '|'.join(s_split[0:(len(s_split)-4)]+list(s[4:8]))
#	outf.write('\t'.join([str(x) for x in s[0:3]+(infos,)])+'\n')