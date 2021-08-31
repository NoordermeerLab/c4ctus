#!/usr/bin/env python


"""
=========================
Steps to generate 4C-Lib
=========================
scripts copied from bbcflib.createlib

"""


import os, random, string, re, sys


def unique_filename_in(path=None):
    """Return a random filename unique in the given path.
    The filename returned is twenty alphanumeric characters which are
    not already serving as a filename in *path*.  If *path* is
    omitted, it defaults to the current working directory.
    """
    if path == None: path = os.getcwd()
    def random_string():
        return "".join([random.choice(string.letters + string.digits) for x in range(20)])
    while True:
        filename = random_string()
        files = [f for f in os.listdir(path) if f.startswith(filename)]
        if files == []: break
    return filename



def parse_fragFile(fragfile,chrom_dict={}):
    """
    Parse fragment file to create segment info bed file and fragment bed file
    """
    segInfoBedFile = "seqInfoBedFile_" + unique_filename_in()
    fragmentBedFile = "fragmentBedFile_" + unique_filename_in()
    o = open(segInfoBedFile,'w')
    obed = open(fragmentBedFile,'w')
    with open(fragfile,'r') as f:
        s = f.next()
        for s in f:
            if re.search('FragIsNotValid',s): continue
            s = s.strip().split('\t')
            chrom = chrom_dict.get(s[1],s[1])
            fragmentInfo = '|'.join(['',s[0],chrom+':'+str(int(s[2])+1)+'-'+s[3],
                                     'indexOfSecondRestSiteOcc='+s[10],
                                     'status='+s[-1],'length='+str(int(s[3])-int(s[2])),
                                     '0','0','0','0'])
            o.write('\t'.join([chrom,s[5],s[6],'type=startSegment'+fragmentInfo])+'\n')
            o.write('\t'.join([chrom,s[8],s[9],'type=endSegment'+fragmentInfo])+'\n')
            row = [chrom,s[2],s[3],'frag'+s[0]]
            obed.write('\t'.join(row)+'\n')
            row[1:3] = [s[5],s[6]]
            obed.write('\t'.join(row)+'_startSeq\n')
            row[1:3] = [s[8],s[9]]
            obed.write('\t'.join(row)+'_endSeq\n')
    o.close()
    obed.close()
    return([segInfoBedFile,fragmentBedFile])




#print "begin"
CD = {}

[s,f] = parse_fragFile(sys.argv[1], CD)

print s
print f
#print "*** done ***"
