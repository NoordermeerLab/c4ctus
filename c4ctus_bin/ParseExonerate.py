#!/usr/bin/env python

"""
===========================
Demultiplexing
Adapted by Benoit Moindrot from: https://github.com/bbcf/bbcflib/blob/master/bbcflib/demultiplex.py
Credit goes to bbcflib
(This script takes the Exonerate output and do the split and file writing)

(Primer file cannot be multiplexed anymore:
To see what multiplexed Primer files were, visit https://github.com/bbcf/bbcflib/blob/master/doc/tutorial_demultiplexing.rst#multiplexed-primer-file)
===========================
"""



import os, string, sys, datetime
import argparse




def split_exonerate(filename,group,minScore,l=30,n=1,trim=True):
    
    #INITIALIZATION
    correction = {}
    files = {}
    filenames = {}
    alignments = {}
    prev_idLine = ''
    
    # Is ambiguous a read with multiple equally good exonerate alignments (i.e. with a same score >= minScore)
    # Will often be empty as it contains only good alignments (>minScore).
    # Other ambiguous but "bad" alignments will be in "unaligned"
    alignments["ambiguous"] = group + "_Myambiguous_"
    files["ambiguous"] = open(alignments["ambiguous"],"w")
    alignments["ambiguous_fastq"] = group + "_Myambiguous-FastQ_"
    files["ambiguous_fastq"] = open(alignments["ambiguous_fastq"],"w")
    # Is unaligned a read with score < minScore (that is score between my_minScore and minScore)
    # my_minScore is the score used for exonerate and is ~ equal to half the length of the primer sequence
    alignments["unaligned"] = group + "_Myunaligned_"
    files["unaligned"] = open(alignments["unaligned"],"w")
    # discarded: if score ok but len(seq) < l (that is the remaining sequence is too short)
    alignments["discarded"] = group + "_Mydiscarded_"
    files["discarded"] = open(alignments["discarded"],"w")
    line_buffer = []
    
    
    def _process_line(line_buffer):
        if len(line_buffer) == 1:
            files[line_buffer[0][3]].write(line_buffer[0][1])
        elif len(line_buffer)>1:
            for buf in sorted(line_buffer):
                files["ambiguous"].write(" ".join(buf[2])+"\n")
            # add the corresponding ambiguous read to the ambiguous_fastq file (only once)
            #files["ambiguous_fastq"].write(" ".join(buf[1])+"\n")
            files["ambiguous_fastq"].write(buf[1])

    with open(filename,"r") as f:
        for s in f:
            if not s[:7] == 'vulgar:': continue
            s = s.strip().split(' ')
            s_split = s[5].split('|')
            key = s_split[0]
            info = s[1].split('_')
            idLine = info[0]
            if idLine != prev_idLine:
                _process_line(line_buffer)
                line_buffer = []
            prev_idLine = idLine
            
            closestBarcode = key
            
            if key not in correction:
                if len(s_split) > 3:
                    correction[key] = len(s_split[1])-len(s_split[3])+n-1
                else:
                    correction[key] = len(s_split[1])+n-1
                
                filenames[key] = group + "_" + key
                files[key] = open(filenames[key],"w")
            k = int(s[3])-int(s[7])+correction[key]
            
            l_linename = len(info[0])
            l_seq = len(info[1])
            full_qual = s[1][int(l_linename)+int(l_seq)+2:int(l_linename)+2*int(l_seq)+2]
            if trim in ['true', '1', 't', 'y', 'yes', 'True', 'TRUE']: ## reads will be trimmed
                seq = info[1][k:l+k] if l+k < l_seq else info[1][k:]; ## if not take the remaining length
                qual = full_qual[k:l+k] if l+k < l_seq else full_qual[k:]
            else: ## reads won't be trimmed
                seq = info[1]
                qual = full_qual
            if int(s[9]) < minScore:
                files["unaligned"].write(" ".join(s)+"\n")
            elif len(seq) >= l/2:
                if key == closestBarcode or closestBarcode == "amb": # if not => not the primer_barcode to which the read comes from + add if ambiguous barcode
                    #line_buffer.append((s[9],"@"+info[0]+"\n"+seq+"***"+info[1]+"\n+\n"+qual+"\n",s,key))
                    line_buffer.append((s[9],"@"+info[0]+"\n"+seq+"\n+\n"+qual+"\n",s,key))
            else: ## will be discarded if remaining read length is too short. Arbitrarly set to < l/2. It used to be < l but will not happen as we now accept l close to the read length.
                files["discarded"].write("@"+info[0]+"\n"+seq+"\n+\n"+qual+"\n")

    _process_line(line_buffer)
    for f in files.itervalues():
        f.close()
    return (filenames,alignments)



if __name__ == "__main__":
    #parse les arguments
    parser = argparse.ArgumentParser(description="Parse exonerate output and split reads based on barcodes")
    parser.add_argument("-f", "--file", required=True, type=str, help="Exonerate output file")
    parser.add_argument("-g", "--group", required=True, type=str, help="Define group [sample / cell line]")
    parser.add_argument("-s", "--minScore", required=True, type=int, help="minScore for parsing exonerate output")
    parser.add_argument("-l", "--length", required=True, type=int, help="Length of the reads to align")
    parser.add_argument("-n", "--FromBase", required=True, type=int, help="Search the primer from base n")
    parser.add_argument("-t", "--trim", required=True, type=str, help="Trim read or not")
    args = vars(parser.parse_args())
    #print args
    
    #print "(" + str(datetime.datetime.now()) + ")    => NOW spliting reads based on Exonerate output"
    #(resSplitExonerate,alignments)=split_exonerate(sys.argv[1],minScore=8,l=30,n=1,trim='False')
    #(resSplitExonerate,alignments)=split_exonerate(sys.argv[1],minScore=int(sys.argv[2]),l=int(sys.argv[3]),n=int(sys.argv[4]),trim=sys.argv[5])
    (resSplitExonerate,alignments)=split_exonerate(args["file"], group=args["group"], minScore=args["minScore"],l=args["length"],n=args["FromBase"], trim=args["trim"])
    #print resSplitExonerate
    for key, value in resSplitExonerate.iteritems():
        print "key={}".format(key) + "\t" + "file={}".format(value)
    for key, value in alignments.iteritems():
        print "key={}".format(key) + "\t" + "file={}".format(value)
    #print alignments
    #print "(" + str(datetime.datetime.now()) + ")    => DONE"