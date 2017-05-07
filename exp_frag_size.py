#!/usr/bin/python3

from bisect import bisect_left
import csv
import re
from tagdigger_fun import reverseComplement, adapters
import argparse
import os

# This script requires a SAM file listing alignments of unique RAD tags to a reference genome.
# Also works on SAM file from UNEAK FASTA file, where there are two tags per marker, named
# 'query' and 'hit'.
# It also requires a FASTA file of the reference genome.
# It outputs a table of expected fragment sizes (without adapters) for each tag.
# These can be used to better understand how fragment size may be biasing read depth.

defcs = 'CTGCAG,CCGG' # default cut sites
defenz = 'PstI-MspI'  # default enzyme pair

maxfragsize = 3000 # trim sequence size for search to save memory (size in bp)

# arguments to the script
parser = argparse.ArgumentParser(description = "TagDigger script for estimating DNA fragment sizes")
parser.add_argument('-s', '--samfile', help = 'SAM file of tags to evaluate', required = True)
parser.add_argument('-g', '--genomefile', help = 'FASTA file of reference genome')
parser.add_argument('-d', '--genome_dir', help = 'Directory with multiple FASTA files of reference genome')
parser.add_argument('-o', '--outfile', help = 'CSV output file', default = 'out.csv')
parser.add_argument('-c', '--cutsites', help = 'Comma-delimited list of restriction sites', default = defcs)
parser.add_argument('-e', '--enzymes', help = 'Name of enzyme pair', default = defenz )
parser.add_argument('-w', '--working_dir', help = 'Directory for reading and writing files')

args = parser.parse_args()
samfile = args.samfile
outputfile = args.outfile

if args.working_dir != None:
    os.chdir(args.working_dir) # change working directories

if (args.genomefile == None and args.genome_dir == None) or \
   (args.genomefile != None and args.genome_dir != None):
    raise Exception("Must provide either one file for reference genome (-g) or directory with multiple files (-d).")

if args.genomefile != None:
    genomefiles = [args.genomefile] # just one genome file
    gfshort = []
else:
    gfshort = os.listdir(args.genome_dir) # genome file names w/o directory
    # genome file names with directory
    genomefiles = [os.path.join(args.genome_dir, x) for x in gfshort]
    gfshort = [x.split('.')[0] for x in gfshort] # remove file extension

# get list of cut sites
if args.enzymes == defenz: # if enzyme pair was left at default, use provided cutsites
    cutsites = [x.strip().upper() for x in args.cutsites.split(',')]
    if not set("".join(cutsites)) <= set('ACGT'):
        raise Exception("Non-ACGT cutsites listed.")
else:
    matchingenz = [x for x in adapters.keys() if x.startswith(args.enzymes)]
    if len(matchingenz) == 0:
        raise Exception("Enzymes {} not found.  See 'adapters' in tagdigger_fun.py.".format(args.enzymes))
    cutsites = [x[0].replace("^","") for x in adapters[matchingenz[0]]]
    if args.cutsites != defcs:
        othercutsites = [x.strip().upper() for x in args.cutsites.split(',')]
        if sorted(cutsites) != sorted(othercutsites):
            raise Exception("Cutsites and enzymes don't match.  Only one of these arguments is needed")

# Begin program

markernames = []   # marker names from first column of SAM file
sequencenames = [] # sequence names from third column of SAM file (e.g. chromosome name)
positions = []     # position of start of marker within its respective sequence
strand = []        # False if reverse complement, True if on forward strand
aligned = []       # True or False for if marker is aligned
tagsizes = []      # length of each tag

# read SAM file
with open(samfile, 'r') as mycon:
    for line in mycon:
        if line[0] == '@': # skip header lines
            continue
        fields = line.split("\t") # split lines into fields
        isQ = fields[0][-8:-3] == "query"
        isH = fields[0][-6:-3] == "hit"
        mrkr = fields[0] # marker or tag name
        seqnm = fields[2] # name of sequence or chromosome
        pos = int(fields[3]) # position of alignment
        flag = bin(int(fields[1]))[2:] # string of bits for flag
        algn = len(flag) < 3 or flag[-3] == '0' # True if strand is aligned
        strd = len(flag) < 5 or flag[-5] == '0' # True if forward strand
        cigar = fields[5] # CIGAR string
        tagsize = len(fields[9])
        if "D" in cigar:
            deletions = sum([int(x) for x in re.findall(r"(\d+)D", cigar)])
        else:
            deletions = 0
        if "I" in cigar:
            insertions = sum([int(x) for x in re.findall(r"(\d+)I", cigar)])
        else:
            insertions = 0
        if "S" in cigar:
            leftpadding = sum([int(x) for x in re.findall(r"(\d+)S", cigar)])
        else:
            leftpadding = 0
        if not strd: # adjust for reverse complement by adding tag size
            pos += tagsize - 1 - insertions + deletions - leftpadding
        else:
            pos -= leftpadding
        if isQ:
            Qmrkr = mrkr.split("_")[0]
            Qseqnm = seqnm
            Qpos = pos
            Qaligned = algn
            Qstrand = strd
        if isH:
            Hmrkr = mrkr.split("_")[0]
            assert Hmrkr == Qmrkr, "UNEAK marker names don't match."
            markernames.append(Hmrkr)
            tagsizes.append(tagsize)
            if algn and seqnm == Qseqnm and pos == Qpos and strd == Qstrand:
                sequencenames.append(seqnm)
                positions.append(pos)
                aligned.append(algn) 
                strand.append(strd)
            else:
                sequencenames.append("*")
                positions.append(0)
                aligned.append(False)
                strand.append(True)
        if (not isQ) and (not isH): # for non-UNEAK formats    
            markernames.append(mrkr)
            tagsizes.append(tagsize)
            sequencenames.append(seqnm)
            positions.append(pos)
            aligned.append(algn) 
            strand.append(strd) 

nMarkers = len(markernames)

# get an order based on the sequence names
listorder = [y for (x,y) in sorted(zip(sequencenames, range(nMarkers)))]
  # -> for each position in the sorted sequences, give the index in the original
seqsort = sorted(sequencenames)

fragsize = ["NA" for i in range(nMarkers)] # for storing the fragment size
outseq = ["" for i in range(nMarkers)] # for the sequence of this fragment
GCcontent = ["NA" for i in range(nMarkers)] # for storing GC content of fragment

# begin reading the FASTA of the genome file
currseqnm = ""
newseqnm = ""
sequence = ""
cnt = 0 # just for counting how many markers found

for i in range(len(genomefiles)):
    with open(genomefiles[i], 'r') as mycon:
        for line in mycon:
            if line[0] == ">": # sequence names
                currseqnm = newseqnm # shift the new one to being the current one
                newseqnm = line[1:].strip() # marker name for the sequence that is to follow
                # check to see if this is the right sequence name or if it should be changed to the file name
                bTest = bisect_left(seqsort, newseqnm)
                if (bTest >= nMarkers or seqsort[bTest] != newseqnm) and len(gfshort) > 0 and \
                   seqsort[bisect_left(seqsort, gfshort[i])] == gfshort[i]: # matches SAM?
                    newseqnm = gfshort[i]
                seqlen = len(sequence)
                if seqlen == 0:
                    continue # skip to next line if there is no sequence yet
                    # search for this sequence name in the alignment info
                b = bisect_left(seqsort, currseqnm) # (designed for large num. scaffolds)
                # loop through any matching markers (or none)
                while b < nMarkers and seqsort[b] == currseqnm: 
                    thisindex = listorder[b] # index of marker in the non-sorted lists
                    if strand[thisindex]:
                        # sequence to search if forward strand
                        subseq = sequence[positions[thisindex]-1:positions[thisindex]+maxfragsize]
                    else:
                        # sequence to search if reverse strand
                        subseq = reverseComplement(sequence[max([0, \
                                  positions[thisindex]-maxfragsize]):positions[thisindex]])
                    size = "NA"
                    for cs in cutsites:
                        # find cut site in the sequence (can't be in tag except at very end)
                        thissize = subseq.find(cs, tagsizes[thisindex] - len(cs)) + len(cs)  
                        if thissize > len(cs) - 1 and (size == "NA" or thissize < size):
                            size = thissize
                    fragsize[thisindex] = size
                    if size != "NA":
                        outseq[thisindex] = subseq[:size]
                        GCcontent[thisindex] = (outseq[thisindex].count('G') + \
                                                outseq[thisindex].count('C'))/size
                    b += 1
                    cnt += 1
                    if cnt % 1000 == 0: # print progress
                        print(cnt)
                # reset sequence
                sequence = ""
            else:
                sequence = sequence + line.strip().upper()
            
# write output file
with open(outputfile, 'w', newline = '') as mycon:
    mywriter = csv.writer(mycon)
    mywriter.writerow(["Marker name", "Sequence name", "Position", "Strand", "Fragment size",
        "Fragment GC content", "Fragment sequence"])
    for i in range(nMarkers):
        mywriter.writerow([markernames[i], sequencenames[i], positions[i], 
                           "forward" if strand[i] else "reverse", fragsize[i], GCcontent[i], outseq[i]])

