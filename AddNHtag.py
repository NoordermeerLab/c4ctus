#!/usr/bin/env python

"""
===========================
Add_NH_tag and generate a BAM file from Bowtie2 output
Benoit   V0.1
[need pysam library - for now, using python on Benoit's MAC
the pysamlibrary is installed on the server -> python]
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


def read_sets(reads,keep_unmapped=False):
    """Groups the alignments in a BAM file by read.

    *reads* should be an iterator over reads, such as the object
     returned by pysam.Samfile.  The SAM/BAM file must be sorted by
     read.  ``read_sets`` removes all unmapped reads, and returns an
     iterator over lists of all AlignedRead objects consisting of the
     same read.
    """
    last_read = None
    for r in reads:
        if (not keep_unmapped) and (r.tid == -1 or r.is_unmapped):
            pass
        elif r.qname != last_read:
            if last_read != None:
                yield accum
            accum = [r]
            last_read = r.qname
        else:
            accum.append(r)
    if last_read != None:
        # We have to check, since if samfile
        # has no alignments, accum is never defined.
        yield accum


def add_nh_flag(samfile, out=None):
    """Adds NH (Number of Hits) flag to each read alignment in *samfile*.

    Scans a BAM file ordered by read name, counts the number of
    alternative alignments reported and writes them to a BAM file
    with the NH tag added.

    If *out* is ``None``, a random name is used.
    """
    if out == None:
        out = unique_filename_in()
    infile = pysam.Samfile(samfile)
    outfile = pysam.Samfile(out, "wb", template=infile)
    for readset in read_sets(infile,keep_unmapped=True):
        nh = sum((r.qlen for r in readset))/(1.0*readset[0].rlen)
        if readset[0].is_paired: nh /= 2
        nh = max(1,int(0.5+nh))
        for read in readset:
            if read.is_unmapped: nh = 0
            read.tags = read.tags+[("NH",nh)]
            outfile.write(read)
    infile.close()
    outfile.close()
    return out




print "Processing " + sys.argv[1]
base    = os.path.basename(sys.argv[1])
OutFile = os.path.splitext(base)[0] + "_NHflaged.bam"
bam = add_nh_flag( sys.argv[1], out = OutFile )
print bam
print "done"