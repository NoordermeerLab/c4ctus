#!/usr/bin/env python

"""
===========================
Restrict the bam files to reads having a NH tag < maxhits
Benoit   V0.1
[need pysam library - for now, using python2 on Benoit's MAC
the pysam library is installed on the server => so python is fine]
===========================
"""



import pysam,os,random,string,sys


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






sorted_bam = sys.argv[1]
maxhits    = int(sys.argv[2])

#reduced_bam = unique_filename_in()
reduced_bam = sys.argv[3]


return_dict = {"fullbam": sorted_bam}

infile = pysam.Samfile( sorted_bam, "rb" )
header = infile.header
#for h in header["SQ"]:
	#if h["SN"] in chromosomes:
	#	h["SN"] = chromosomes[h["SN"]]["name"]
outfile = pysam.Samfile( reduced_bam, "wb", header=header )
for read in infile:
	nh = dict(read.tags).get('NH',1)
	if nh < 1:
		continue
	if nh <= maxhits:
		outfile.write(read)
		#print "nh = " + str(nh) + "   =>   read is kept"
outfile.close()
infile.close()
#        reduced_bam = sort_bam.nonblocking(ex, bam2, via=via).wait()
#index2 = index_bam(ex, reduced_bam)
#        reduced_bam = add_and_index_bam( ex, bam2, set_file_descr(name+"filtered.bam",**bam_descr) )
return_dict['filteredbam'] = reduced_bam
#return_dict['stats'] = full_stats
#return return_dict

print return_dict