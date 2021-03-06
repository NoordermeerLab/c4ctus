#!/usr/bin/env python

# Download from: https://raw.githubusercontent.com/bbcf/bbcfutils/master/Python/fastqToFasta.py
# 
# Adapted to not require bbcflib.common by Benoit Moindrot




import sys, optparse, pysam, os, random, string



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



usage = sys.argv[0]+" [OPTIONS]"
description = "Conversion from fastq to fasta."
opts = (("-i", "--input", "Fastq filename.", {'default':None}),
        ("-o", "--output", "Output fasta filename (default Random name).", 
         {'default':None}),
        ("-n", "--start", "Base number to start from", {'default': 1}),
        ("-x", "--length", "Length of reads in output", {'default': 22}),
        ("-d", "--debug", "Print debugging info", 
         {'action':"store_true",'default': False}))



class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg


def main():
    parser = None
    try:
        parser = optparse.OptionParser(usage=usage, description=description)
        for opt in opts:
            if len(opt) == 4:
                parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
            elif len(opt) == 3:
                parser.add_option(opt[0],help=opt[1],**opt[2])
        (opt, args) = parser.parse_args()

        if not(opt.input and os.path.exists(opt.input)):
            raise Usage("Please provide a fastq file")

        if opt.debug:
            print("""
fastqToFasta.py
i=%s
n=%i
x=%i
""" %(opt.input,n,x))

        fq = pysam.FastqFile(opt.input)
        faFile = opt.output or unique_filename_in()
        rlen = int(opt.length)
        rskip = int(opt.start)-1
        fa = open(faFile,"w")    
        for i,s in enumerate(fq):
            seq = s.sequence[rskip:(rskip+rlen)]
            header = "_".join([s.name,s.sequence,s.quality])
            fa.write(">"+header+"\n"+seq+"\n")
        fq.close()
        fa.close()
        print "Converting fastq to fasta= DONE  (file = " + faFile + ")"

    except Usage, err:
        print >>sys.stderr, '\n',err.msg,'\n'
        if parser: parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main())

