# Functions and restriction enzyme data for the Python program TagDigger,
# created by Lindsay V. Clark.

import gzip
import csv
import hashlib
import os
import sys
import bisect
import re

# confirm that Python 3 is being used.
if sys.version_info.major < 3 or (sys.version_info.major == 3 and sys.version_info.minor < 3):
    print("Python 3.3 or higher required.  You appear to be using an older version of Python.")
    raise Exception("Python 3 required.")

# dictionary of restriction enzyme cut sites, including only what will show
# up after the barcode.
enzymes = {'ApeKI': 'CWGC', 'EcoT22I': 'TGCAT', 'NcoI': 'CATGG',
           'NsiI': 'TGCAT', 'PstI': 'TGCAG', 'SbfI': 'TGCAGG', 'None': ''}

# dictionary of possible adapter sequences to find.
# Restriction site listed first, with ^ indicating the end of genomic
# sequence that we would expect to find.  Then top strand adapter sequence,
# not including the restriction site overhang.  [barcode] indicates where the
# reverse complement of the barcode should be located.
adapters = {'PstI-MspI-Hall': [('CCG^G', # Sacks lab, designed by Megan Hall (I think)
                                'CTCAGGCATCACTCGATTCCTCCGTCGTATGCCGTCTTCTGCTTG'),
                               ('CTGCA^G',
                                '[barcode]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')],
            'NsiI-MspI-Hall': [('CCG^G', # Sacks lab, NsiI, with above adapters
                                'CTCAGGCATCACTCGATTCCTCCGTCGTATGCCGTCTTCTGCTTG'),
                               ('ATGCA^T',
                                '[barcode]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')],
            'PstI-MspI-Clark': [('CCG^G', # Sacks lab, altered for full P7 sequence
                                'CTCAGGCATCACTCGATTCCTATCTCGTATGCCGTCTTCTGCTTG'),
                               ('CTGCA^G',
                                '[barcode]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')],
            'NsiI-MspI-Clark': [('CCG^G', # Sacks lab, NsiI, with above adapters for full P7 sequence
                                'CTCAGGCATCACTCGATTCCTATCTCGTATGCCGTCTTCTGCTTG'),
                               ('ATGCA^T',
                                '[barcode]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')],
            'PstI-MspI-Poland': [('CCG^G', # DOI: 10.1371/journal.pone.0032253
                                  'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG'),
                                 ('CTGCA^G',
                                  '[barcode]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')]
            }

# codes for ambiguous nucleotides
IUPAC_codes = {frozenset('AG'): 'R', frozenset('CT'): 'Y',
               frozenset('GT'): 'K', frozenset('AC'): 'M',
               frozenset('CG'): 'S', frozenset('AT'): 'W',
               frozenset('CGT'): 'B', frozenset('AGT'): 'D',
               frozenset('ACT'): 'H', frozenset('ACG'): 'V',
               frozenset('ACGT'): 'N',
               frozenset('A'): 'A', frozenset('C'): 'C',
               frozenset('G'): 'G', frozenset('T'): 'T'}

# Function definitions.
def combine_barcode_and_cutsite(barcodes, cutsite):
    '''Add restriction cut sites to the end of each barcode in a list,
    to make searching for barcode + cutsite more efficient.  Convert
    to upper case.'''
    assert all([set(barcode.upper()) <= set('ACGT') for
                barcode in barcodes]), "Non-ACGT barcode."
    assert set(cutsite.upper()) <= set('ACGT'), "Invalid cut site."
    # Note that cut sites have already been enumerated and therefore should not
    # have ambiguous letters at this point.
    return [(barcode + cutsite).upper() for barcode in barcodes]

def tree_one_level(num_seq):
    '''Take a set of sequences, sort them by the first letter, return
     a set of lists.  num_seq is a list of lists, where each list has
     the sequence and its index number.'''
    # if we have gotten to the end of a sequence
    if len(num_seq[0][0]) == 0:
        return num_seq[0]
    # if we need to split sequences up by the first letter
    else:
        mysubtree = [['A',[]],['C',[]],['G',[]],['T',[]]]
        for s in num_seq:
            assert len(s[0]) > 0, "Problematic sequence: {}.  Likely due to overlapping tags.".format(s[1])
            myind = "ACGT".find(s[0][0])
            s[0] = s[0][1:]
            mysubtree[myind][1].append(s)
        return mysubtree

def tree_recursive(num_seq):
    '''Build tree recursively for indexing sequences (either barcodes or
    tags) rapidly.'''
    mysubtree = tree_one_level(num_seq)
    if len(mysubtree[0]) > 0:
        for i in range(4):
            if len(mysubtree[i][1]) > 0:
                mysubtree[i][1] = tree_recursive(mysubtree[i][1])
    return mysubtree

def build_sequence_tree(sequences, numseq):
    '''Give each sequence a number to use for indexing, then build tree.
    numseq is maximum number to use for indexing, which is len(sequences)
    unless there are multiple cutsites (like with ApeKI).'''
    myindex = 0
    seqlist = []
    for s in sequences:
        seqlist.append([s, myindex])
        myindex += 1
        if myindex == numseq:
            myindex = 0
    if(numseq == 1 and sequences == ['']):
        result = [['A', ['', 0]], ['C', ['', 0]], ['G', ['', 0]], ['T', ['', 0]]]
    else:
        result = tree_recursive(seqlist)
    return(result)

def sequence_index_lookup(sequence, seqtree):
    '''Lookup a sequence in an index tree, and get the index.  Return -1 if
       the sequence is not in the index tree.'''
    # If we are out of sequence
    if len(sequence) == 0:
        return -1
    # If we get an N or something else not in ACGT:
    if sequence[0] not in {'A', 'C', 'G', 'T'}:
        return -1
    # If we do have a normal nucleotide, get its index for the tree
    myind = "ACGT".find(sequence[0])
    # If we have gotten to the point where we can tell the sequence
    # does not exist in this set:
    if len(seqtree[myind][1]) == 0:
        return -1
    # If we have gotten to the end of the sequence and can retrieve the index:
    if len(seqtree[myind][1][0]) == 0:
        return seqtree[myind][1][1]
    # otherwise go to the next nucleotide
    return sequence_index_lookup(sequence[1:], seqtree[myind][1])

def enumerate_cut_sites(cutsite):
    '''Given an ambiguous cut site written in IUPAC codes, generate all
    possible cut sites.'''
    cutsites = [cutsite]
    while cutsites[0].find('R') > -1:
        newcutsites = [x.replace('R', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('R', 'G', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('Y') > -1:
        newcutsites = [x.replace('Y', 'C', 1) for x in cutsites]
        newcutsites += [x.replace('Y', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('K') > -1:
        newcutsites = [x.replace('K', 'G', 1) for x in cutsites]
        newcutsites += [x.replace('K', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('M') > -1:
        newcutsites = [x.replace('M', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('M', 'C', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('S') > -1:
        newcutsites = [x.replace('S', 'C', 1) for x in cutsites]
        newcutsites += [x.replace('S', 'G', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('W') > -1:
        newcutsites = [x.replace('W', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('W', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('B') > -1:
        newcutsites = [x.replace('B', 'C', 1) for x in cutsites]
        newcutsites += [x.replace('B', 'G', 1) for x in cutsites]
        newcutsites += [x.replace('B', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('D') > -1:
        newcutsites = [x.replace('D', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('D', 'G', 1) for x in cutsites]
        newcutsites += [x.replace('D', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('H') > -1:
        newcutsites = [x.replace('H', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('H', 'C', 1) for x in cutsites]
        newcutsites += [x.replace('H', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('V') > -1:
        newcutsites = [x.replace('V', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('V', 'C', 1) for x in cutsites]
        newcutsites += [x.replace('V', 'G', 1) for x in cutsites]
        cutsites = newcutsites
    while cutsites[0].find('N') > -1:
        newcutsites = [x.replace('N', 'A', 1) for x in cutsites]
        newcutsites += [x.replace('N', 'C', 1) for x in cutsites]
        newcutsites += [x.replace('N', 'G', 1) for x in cutsites]
        newcutsites += [x.replace('N', 'T', 1) for x in cutsites]
        cutsites = newcutsites
    return cutsites

def find_tags_fastq(fqfile, barcodes, tags, cutsite="TGCAG",
                    maxreads=5000000000, tassel_tagcount = False):
    '''Make indexing trees for barcodes, cut sites, and tags, then go
    through a FASTQ file and count up instances of barcode*tag combos.'''

    # asserts
    assert all([set(barcode.upper()) <= set('ACGT') for
                barcode in barcodes]), "Non-ACGT barcode."
    cutsite = cutsite.upper()
    assert set(cutsite) <= set('ACGTNRYKMSWBDHV'), "Invalid cut site."
    tags = [tag.upper() for tag in tags]
    assert all([set(tag) <= set('ACGT') for tag in tags]), \
           "Non-ACGT tag."

    # Get cut site length
    cutlen = len(cutsite)
    # Get lengths of all barcodes + cut site
    barcutlen = [len(x) + cutlen for x in barcodes]
    # Get number of barcodes
    barnum = len(barcodes)
    # Make multiple cut sites if there are any ambiguous bases.
    cutsites = enumerate_cut_sites(cutsite)
    # Combine barcodes with cut sites.
    barcut = []
    for cut in cutsites:
        barcut += combine_barcode_and_cutsite(barcodes, cut)
    # Make indexing tree for barcodes + cut sites.
    barcuttree = build_sequence_tree(barcut, barnum)

    # If we see every tag starting with a cut site:
    if set([x[:cutlen] for x in tags]).issubset(set(cutsites)):
        # If there is only one possible cut site, remove it from tags since
        # we already checked for it when checking for barcodes.
        if(len(cutsites) == 1):
            tags = [x[cutlen:] for x in tags]
        else:
        # If there are multiple cutsites, we need to keep them as part of the
        # tags (they may contain a relevant SNP) and adjust the starting
        # position for where we will read tags in the FASTQ sequences.
            barcutlen = [x - cutlen for x in barcutlen]
    # Make indexing tree for tags.
    tagtree = build_sequence_tree(tags, len(tags))

    # Set up matrix to contain tag counts.  First index is barcode and second
    # index is tag.
    mycounts = [[0 for x in range(len(tags))] for x in range(barnum)]

    # Open connection to FASTQ file.
    if fqfile[-2:].lower() == 'gz':
        fqcon = gzip.open(fqfile, 'rt')
    else:
        fqcon = open(fqfile, 'r')
    # Loop through sequence reads and search for tags.
    try:
        readscount = 0
        barcutcount = 0
        tagcount = 0
        lineindex = 0
        for line in fqcon:
            if lineindex % 4 == 0 and tassel_tagcount:
                # extract counts if the file is converted from a TASSEL tagCount file
                thesecounts = int(line[line.find("count=")+6:].strip())
            if lineindex % 4 == 1: # lines with sequence
                readscount += 1
                line1 = line.strip().upper()
                barindex = sequence_index_lookup(line1, barcuttree)
                if barindex > -1: # if barcode and cut site have a match
                    barcutcount += 1
                    tagindex = sequence_index_lookup(line1[barcutlen[barindex]:],
                                                     tagtree)
                    if tagindex > -1: # if the sequence matches an expected tag
                        tagcount += 1
                        if tassel_tagcount:
                            mycounts[barindex][tagindex] += thesecounts
                        else:
                            mycounts[barindex][tagindex] += 1
                if readscount % 1000000 == 0:
                    print(fqfile)
                if readscount % 50000 == 0:
                    print("Reads: {0} With barcode and cut site: {1} With tag: {2}".format(readscount, barcutcount, tagcount))
                if readscount >= maxreads: # stop loop if we get to max reads
                    break
            lineindex += 1
    finally:
        fqcon.close()
    return mycounts

def isFastq(filename):
    '''Check the first four lines of a file, return 1 if is it an unzipped
       FASTQ file, 2 if it is gzipped FASTQ file, and 0 if it is neither or
       if it can't be opened.'''
    try:
        # determine whether to open as zipped or not
        if filename[-2:].lower() == 'gz':
            mycon = gzip.open(filename, 'rt')
            myresult = 2
        else:
            mycon = open(filename, 'r')
            myresult = 1
    except IOError:
        myresult = 0
    else:
        try:
            thisline = mycon.readline() # first line is comment
            if thisline[0] != '@':
                myresult = 0
            thisline = mycon.readline() # second line is sequence
            if not set(thisline.strip()) <= {'A', 'C', 'G', 'T', 'N',
                                             'a', 'c', 'g', 't', 'n'}:
                myresult = 0
            thisline = mycon.readline() # third line is comment
            if thisline[0] != '+':
                myresult = 0
        finally:
            mycon.close()
    return myresult

def readBarcodeKeyfile(filename, forSplitter=False):
    '''Read in a csv file containing file names, barcodes, and sample names,
       and output a dictionary containing this information.
       forSplitter: indicates whether this is for the barcode splitter, as
       opposed to the tag counter.'''
    try:
        with open(filename, 'r', newline='') as mycon:
            mycsv = csv.reader(mycon)
            rowcount = 0
            result = dict()
            for row in mycsv:
                if rowcount == 0:
                    # read header row
                    if forSplitter:
                        fi = row.index("Input File")
                        bi = row.index("Barcode")
                        si = row.index("Output File")
                    else:
                        fi = row.index("File")
                        bi = row.index("Barcode")
                        si = row.index("Sample")
                else:
                    f = row[fi].strip() # file name
                    b = row[bi].strip().upper() # barcode
                    s = row[si].strip() # sample
                    # skip blank lines
                    if f == "" and b == "" and s == "":
                        continue
                    # check row contents
                    if f == "":
                        raise Exception("Blank cell found where file name should be in row {}.".format(rowcount + 1))
                    if s == "":
                        raise Exception("Blank cell found where sample name should be in row {}.".format(rowcount + 1))
                    if not set(b) <= {'A', 'C', 'G', 'T'}:
                        raise Exception("{0} in row {1} is not a valid barcode.".format(b, rowcount + 1))
                    # check whether this file is already in dict, add if necessary
                    if f not in result.keys():
                        result[f] = [[], []] # first list for barcodes, second for samples
                    if b in result[f][0]:
                        raise Exception("Each barcode can only be present once for each file.")
                    # add the barcode and sample
                    result[f][0].append(b)
                    result[f][1].append(s)
                rowcount += 1
        # check that no output file names are duplicated
        if forSplitter:
            alloutfile = set() # all unique output files
            totout = 0 # total number of output files
            for k in result.keys():
                alloutfile.update(result[k][1])
                totout += len(result[k][1])
            if len(alloutfile) < totout:
                raise Exception("All output files must have unique names for barcode splitter.")
    except IOError:
        print("Could not read file {}.".format(filename))
        result = None
    except ValueError:
        if forSplitter:
            print("File header needed containing 'Input File', 'Barcode', and 'Output File'.")
        else:
            print("File header needed containing 'File', 'Barcode', and 'Sample'.")
        result = None
    except Exception as err:
        print(err.args[0])
        result = None
    return result

def compareTags(taglist, trim = True):
    '''Find SNP alleles within a set of sequence tags representing the
       same locus.'''
    assert type(taglist) is list, "taglist must be list."
    assert all([set(t) <= set('ATCG') for t in taglist]), \
           "taglist must be a list of ACGT strings."
    # make sure all tags are same length
    if len(set([len(t) for t in taglist])) > 1:
        if trim: # trim down to minimum length
            minlen = min([len(t) for t in taglist])
            taglist = [tag[:minlen] for tag in taglist]
        else: # pad out to maximum length
            maxlen = max([len(t) for t in taglist])
            taglist = [tag.ljust(maxlen, 'N') for tag in taglist]
    # get index(es) and nucleotides for variable sites
    result =[(i, [t[i] for t in taglist]) for i in range(len(taglist[0])) \
                 if len(set([t[i] for t in taglist if t[i] != 'N'])) > 1]
    return result

def readTags_UNEAK_FASTA(filename, toKeep = None):
    '''Read in pairs of tags from a FASTA file output by the TASSEL UNEAK
       pipeline. filename is FASTA file, and toKeep is a list of tag names
       to retain.'''
    linecount = 0
    namelist = [] # for storing tag names
    seqlist = []  # for storing tag sequences
    try:
        with open(filename, mode='r') as mycon:
            for line in mycon:
                if linecount % 4 == 0: # lines with first tag name
                    if line[:3] != ">TP":
                        raise Exception("Line {0} of {1} does not start with '>TP'.".format(linecount+1, filename))
                    tagname1 = line[1:line.rfind("_")]
                    # get length of real tag (some are padded with A's at the end)
                    taglength1 = int(line[line.rfind("_")+1:].strip())
                if linecount % 4 == 1: # lines with first tag sequence
                    seq1 = line.strip().upper()
                    seq1 = seq1[:taglength1]
                    if not set(seq1) <= set('ACGT'):
                        raise Exception("Line {0} is not ACGT sequence.".format(linecount+1))
                    if seq1 in seqlist:
                        raise Exception("Non-unique sequence found: line {0}.".format(linecount+1))
                if linecount % 4 == 2: # lines with second tag name
                    if line[:3] != ">TP":
                        raise Exception("Line {0} of {1} does not start with '>TP'.".format(linecount+1, filename))
                    tagname2 = line[1:line.rfind("_")]
                    if tagname1[:tagname1.find("_")] != tagname2[:tagname2.find("_")]:
                        raise Exception("Tag name in line {0} does not match tag name in line {1}.".format(linecount+1, linecount - 1))
#                    if taglength != int(line[line.rfind("_")+1:].strip()):
#                        raise Exception("Tag length in line {0} does not match tag length in line {1}.".format(linecount+1, linecount - 1))
                    ## Note: some tag pairs output by UNEAK do have different lengths.  Polymorphism in second cut site?
                    taglength2 = int(line[line.rfind("_")+1:].strip())
                if linecount % 4 == 3: # lines with second tag sequence
                    seq2 = line.strip().upper()
                    seq2 = seq2[:taglength2]
                    if not set(seq2) <= set('ACGT'):
                        raise Exception("Line {0} is not ACGT sequence.".format(linecount+1))
                    if seq2 in seqlist:
                        raise Exception("Non-unique sequence found: line {0}.".format(linecount+1))
                    # determine whether to skip this one if just doing a subset
                    if toKeep != None and tagname1[:tagname1.find("_")] not in toKeep:
                        linecount +=1
                        continue
                    # if tags are different lengths, make sure that they can be distinguished
                    minlen = min(taglength1, taglength2)
                    if taglength1 != taglength2 and seq1[:minlen] == seq2[:minlen]:
                        print("{} skipped because tags cannot be distinguished.".format(tagname1[:tagname1.find("_")]))
                        linecount += 1
                        continue

                    # determine differences between sequences
                    diff = compareTags([seq1, seq2])
                    # add nucleotide to tag name
                    tagname1 += "_"
                    tagname1 += diff[0][1][0]
                    tagname2 += "_"
                    tagname2 += diff[0][1][1]
                    # determine which should be the first and second allele, using
                    # alphabetical order to match hapMap2numeric.
                    if diff[0][1][0] < diff[0][1][1]:
                        tagname1 += "_0"
                        tagname2 += "_1"
                    else:
                        tagname1 += "_1"
                        tagname2 += "_0"
                    # add names and sequences to lists
                    namelist.extend([tagname1, tagname2])
                    seqlist.extend([seq1[:minlen], seq2[:minlen]])
                linecount += 1
        result = [namelist, seqlist]
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    except Exception as err:
#        print("Line {}".format(linecount +1))
        print(err.args[0])
        result = None
    return result

def readTags_Rows(filename, toKeep = None):
    '''Read a CSV of sequence tags, where each row has a marker name, an
       allele name, and a tag sequence.'''
    rowcount = 0
    namelist = []
    seqlist = []
    try:
        with open(filename, mode='r') as mycon:
            mycsv = csv.reader(mycon)
            for row in mycsv:
                if rowcount == 0: # header row
                    if not {"Marker name", "Allele name", "Tag sequence"} <= set(row):
                        raise Exception("Need 'Marker name', 'Allele name', and 'Tag sequence' in header row.")
                    mi = row.index("Marker name")
                    ai = row.index("Allele name")
                    ti = row.index("Tag sequence")
                else: # rows after header row
                    mname = row[mi].strip() # marker name
                    if '_' in mname:
                        raise Exception("Marker {}: marker names cannot contain underscores.".format(mname))
                    if toKeep != None and mname not in toKeep:
                        rowcount += 1
                        continue # skip this marker if it is not in the list of ones we want.
                    allele = row[ai].strip() # allele name
                    tag = row[ti].upper().strip() # tag sequence
                    if not set(tag) <= set('ACGT'):
                        raise Exception("Tag sequence not formatted as ACGT in row {}.".format(rowcount+1))
                    if tag in seqlist:
                        raise Exception("Non-unique sequence found: line {0}.".format(rowcount+1))
                    namelist.append(mname + '_' + allele)
                    seqlist.append(tag)
                rowcount += 1
        result = [namelist, seqlist]
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    except Exception as err:
        print(err.args[0])
        result = None
    return result

def readTags_Columns(filename, toKeep = None):
    '''Read a CSV of sequence tags, where each row has a marker name and
       two tag sequences representing two alleles.'''
    rowcount = 0
    namelist = []
    seqlist = []
    try:
        with open(filename, mode='r') as mycon:
            mycsv = csv.reader(mycon)
            for row in mycsv:
                if rowcount == 0:
                    if not {"Marker name", "Tag sequence 0", "Tag sequence 1"} <= set(row):
                        raise Exception("Need 'Marker name', 'Tag sequence 0', and 'Tag sequence 1' in header row.")
                    mi = row.index("Marker name")
                    ti0 = row.index("Tag sequence 0")
                    ti1 = row.index("Tag sequence 1")
                else: # all rows after header row
                    mname = row[mi].strip() # marker name
                    if '_' in mname:
                        raise Exception("Marker {}: marker names cannot contain underscores.".format(mname))
                    if toKeep != None and mname not in toKeep:
                        rowcount +=1
                        continue # skip this marker if it is not in the list of ones we want.
                    tag0 = row[ti0].upper().strip() # get tag sequences
                    tag1 = row[ti1].upper().strip()
                    if not set(tag0 + tag1) <= set('ACGT'):
                        raise Exception("Tag sequence not formatted as ACGT in row {}.".format(rowcount+1))
                    if (tag0 in seqlist) or (tag1 in seqlist):
                        raise Exception("Non-unique sequence found: line {0}.".format(rowcount+1))
                    seqlist.append(tag0)
                    seqlist.append(tag1)
                    # make allele names
                    diff = compareTags([tag0, tag1])
                    t0SNP = ''.join([s[1][0] for s in diff])
                    t1SNP = ''.join([s[1][1] for s in diff])
                    namelist.append(mname + '_' + t0SNP + '_0')
                    namelist.append(mname + '_' + t1SNP + '_1')
                rowcount += 1
            result = [namelist, seqlist]
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    except Exception as err:
        print(err.args[0])
        result = None
    return result

def readTags_Merged(filename, toKeep = None, allowDuplicates=False):
    '''Read a CSV of sequence tags, where each row has a marker name and a
       tag sequence with polymorphic nucleotides in square brackets separated
       by a forward slash.  Because this is the format for the SNP database
       produced by the tag manager, and that database might have multiple
       markers that share alleles (since the tag manager currently treats all
       markers as biallelic), this function differs from the other readTags
       functions in how it handles duplicate tags.'''
    rowcount = 0
    namelist = []
    seqlist = []
    try:
        with open(filename, mode='r') as mycon:
            mycsv = csv.reader(mycon)
            for row in mycsv:
                if rowcount == 0:
                    if not {"Marker name", "Tag sequence"} <= set(row):
                        raise Exception("Need 'Marker name' and 'Tag sequence' in header row.")
                    mi = row.index("Marker name")
                    ti = row.index("Tag sequence")
                else:
                    if not set('[/]') < set(row[ti]):
                        raise Exception("Characters '[/]' not found in row {}.".format(rowcount+1))
                    mname = row[mi].strip() # marker name
                    if '_' in mname:
                        raise Exception("Marker {}: marker names cannot contain underscores.".format(row[mi]))
                    if toKeep != None and mname not in toKeep:
                        rowcount += 1
                        continue # skip this marker if it is not in the list of ones we want.
                    # find delimiting characters within sequence
                    p1 = row[ti].find('[')
                    #p2 = row[ti].find('/')
                    p3 = row[ti].find(']')
                    # recreate the tag sequences
                    subtags = row[ti][p1+1:p3].split('/') # list of versions of the variable region
                    subtags = [x.strip().upper() for x in subtags]
                    tags = [(row[ti][:p1] + x + row[ti][p3+1:]).upper().strip().replace('-','') for x in subtags]
                    
                    if not allowDuplicates and any([x in seqlist for x in tags]):
                        print("Non-unique sequence found: line {0}.".format(rowcount+1))
                        print("Marker {} skipped.".format(mname))
                        rowcount += 1
                        continue
                    seqlist.extend(tags)
                    if not all([set(x) <= set('ACGT') for x in tags]):
                        raise Exception("Tag sequence not formatted correctly in row {}.".format(rowcount+1))
                    # generate the allele names
                    namelist.extend(["{}_{}_{}".format(mname, subtags[i], i) for i in range(len(tags))])
                rowcount += 1
        result = [namelist, seqlist]
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    except Exception as err:
        print(err.args[0])
        result = None
    return result

def readTags_Stacks(tagsfile, snpsfile, allelesfile, toKeep = None, binaryOnly=False):
    '''Read tags from the catalog format produced by Stacks.'''
    try:
        alltags = dict() # keys are locus numbers, values are sequences
        with gzip.open(tagsfile, mode = 'rt') if tagsfile.endswith('.gz') else open(tagsfile, mode = 'r') as mycon:
            tr = csv.reader(mycon, delimiter='\t')
            for row in tr:
                if row[0].startswith("#"):
                    continue  # skip comment line
                if toKeep == None or row[2] in toKeep:
                    alltags[row[2]] = row[9]
        alleles = list() # tuples, where first item is locus number and second is haplotype
        with gzip.open(allelesfile, mode = 'rt') if allelesfile.endswith('.gz') else open(allelesfile, mode = 'r') as mycon:
            ar = csv.reader(mycon, delimiter='\t')
            for row in ar:
                if row[0].startswith("#"):
                    continue
                if toKeep == None or row[2] in toKeep:
                    alleles.append((row[2], row[3]))
        positions = dict() # keys are locus numbers, values are lists of variant positions
        with gzip.open(snpsfile, mode = 'rt') if snpsfile.endswith('.gz') else open(snpsfile, mode = 'r') as mycon:
            sr = csv.reader(mycon, delimiter='\t')
            for row in sr:
                if row[0].startswith("#"):
                    continue
                if toKeep == None or row[2] in toKeep:
                    if row[2] in positions.keys():
                        positions[row[2]].append(int(row[3]))
                    else:
                        positions[row[2]] = [int(row[3])]
        namelist = []
        seqlist = []
        # generate tag sequences
        for a in alleles:
            thislocus = a[0]
            consensusseq = alltags[thislocus]
            if len(a[1]) == 0: # non-variable tags
                outseq = consensusseq
            else:
                thesepos = positions[thislocus]
                outseq = consensusseq[:thesepos[0]] # first non-variable portion
                for i in range(len(a[1])):
                    outseq = outseq + a[1][i]
                    if i+1 == len(a[1]):
                        outseq = outseq + consensusseq[thesepos[i]+1:]
                    else:
                        outseq = outseq + consensusseq[thesepos[i]+1:thesepos[i+1]]
            outseq = outseq.upper()
            if not set(outseq) <= set('ACGT'):
                print("{}_{} skipped for having non-ACGT nucleotides.".format(a[0], a[1]))
            else:
                namelist.append(a[0] + '_' + a[1]) # name by locus number and haplotype
                seqlist.append(outseq)
        # cleanup
        if binaryOnly:
            markers = extractMarkers(namelist)
            newnamelist = []
            newseqlist = []
            for i in range(len(markers[0])):
                if len(markers[1][i][0]) != 2:
                    continue
                index1 = markers[1][i][1][0]
                index2 = markers[1][i][1][1]
                # give 0 or 1 designation based on alphabetical order
                if markers[1][i][0][0] < markers[1][i][0][1]:
                    newname1 = namelist[index1] + '_0'
                    newname2 = namelist[index2] + '_1'
                else:
                    newname1 = namelist[index1] + '_1'
                    newname2 = namelist[index2] + '_0'
                newnamelist.extend([newname1, newname2])
                newseqlist.extend([seqlist[index1], seqlist[index2]])
            namelist = newnamelist
            seqlist = newseqlist
        return [namelist, seqlist]

    except IOError:
        print("Files not readable.")
        return None
    except IndexError:
        print("Files in wrong format.")
        return None
    except ValueError:
        print("Files in wrong format.")
        return None
    except KeyError:
        print("Locus names not matching properly.")
        return None
    except Exception as err:
        print(err.args[0])
        return None

def readTags_TASSELSAM(filename, toKeep=None, binaryOnly=False, writeMarkerKey=False,
                       keyfilename = None):
    '''Read tag sequences from a SAM file and match marker names using the
       same conventions as TASSEL-GBSv2.
       toKeep: an optional list of marker names to keep, in TASSEL-GBSv2 format.
       binaryOnly: boolean indicating whether to only keep markers with two tags.
       writeMarkerKey: boolean indicating whether to write a key of TASSEL SNP names vs.
           tagDigger marker names.
       keyfilename: name for marker key file.'''
    assert (not writeMarkerKey) or keyfilename != None, "keyfilename needed."
    namelist = []
    seqlist = []
    tempseq = dict() # keys are marker names, values are lists of tags
    numdig = 0 # number of digits for indicating alignment position
    markerkey = [] # each item will be a tuple containing the TASSEL name for the SNP followed by
                    # the tagdigger name for the marker.
    try:
        with open(filename, mode='r') as mycon:
            for line in mycon:
                if line[0:3] == '@SQ':
                    # length of number indicating chromosome size
                    chrsize = len(line.split()[2][3:])
                    if chrsize > numdig:
                        numdig = chrsize
                    continue
                elif line[0] == '@':
                    continue

                mycolumns = line.split()
                myflags = int(mycolumns[1])
                # skip if no alignment (4 flag)
                if myflags - 4 in {0, 1, 2, 8, 16, 32, 64, 128}:
                    continue

                # chromosome name
                chrom = mycolumns[2] #.upper()
                # eliminate underscores in chromosome names
                chrom = chrom.replace('_', '*')
                # position
                pos = int(mycolumns[3])
                # strand
                if myflags - 16 in {0, 1, 2, 8, 32, 64, 128}:
                    strand = "bot"
                else:
                    strand = "top"
                # sequence
                sequence = mycolumns[9]
                if strand == 'bot':
                    sequence = reverseComplement(sequence)
                    # adjust position to begin at cut site
                    cigar = mycolumns[5]
                    deletions = sum([int(x[:-1]) for x in re.findall('\d+D', cigar)])
                    insertions = sum([int(x[:-1]) for x in re.findall('\d+I', cigar)])
                    pos = pos + len(sequence) - insertions + deletions - 1

                # make marker name and add to dictionary
                marker ="{}-{:0>{width}}-{}".format(chrom, pos, strand, width=numdig)
                if marker in tempseq.keys():
                    # check that an overlapping version of this sequence isn't already there (keep shorter version).
                    # can happen rarely with restriction site polymorphism.
                    addthisseq = True
                    tempseq[marker] = [ts for ts in tempseq[marker] if not ts.startswith(sequence)]
                    for exstseq in tempseq[marker]:
                        if sequence.startswith(exstseq):
                            addthisseq = False # this one is longer, don't add it
                    if addthisseq:
                        tempseq[marker].append(sequence)
                else:
                    tempseq[marker] = [sequence]

        allMarkers = sorted(tempseq.keys())
        for m in allMarkers:
            thesetags = tempseq[m]
            ntags = len(thesetags)
            if binaryOnly and ntags != 2:
                continue
            diff = compareTags(thesetags, trim = False)

            if toKeep != None or writeMarkerKey: # check if marker is in list to be retained
                markerinfo = m.split('-')
                chrom = markerinfo[0].upper()
                if chrom.startswith("CHROMOSOME"):
                    chrom = chrom[10:]
                if chrom.startswith("CHR"):
                    chrom = chrom[3:]
                # put S at the beginning of chromosome number ## is this always right?
                if chrom[0] in set('0123456789'):
                    chrom = 'S' + chrom
                pos = int(markerinfo[1]) # position of leftmost nucleotide
                if markerinfo[2] == 'top':
                    snppos = [pos + d[0] for d in diff]
                else:
                    snppos = [pos - d[0] for d in diff]
                # generate SNP names in TASSEL format
                possiblenames = ['{}_{}'.format(chrom, p) for p in snppos]
                # skip if none of these SNP names were in the list to keep
                if toKeep != None and all([p not in toKeep for p in possiblenames]):
                    continue
                if writeMarkerKey:
                    # add names to list of markers to write.
                    for p in possiblenames:
                        markerkey.append((p, m))
            allelenames = [''.join([d[1][i] for d in diff]) for i in range(ntags)]
            thesetagnames = [m + '_' + a for a in allelenames]
            if binaryOnly and allelenames[0] < allelenames[1]:
                thesetagnames[0] += '_0'
                thesetagnames[1] += '_1'
            if binaryOnly and allelenames[1] < allelenames[0]:
                thesetagnames[1] += '_0'
                thesetagnames[0] += '_1'
            namelist.extend(thesetagnames)
            seqlist.extend(thesetags)
        if len(namelist) == 0:
            raise Exception("No markers output; is list of markers to keep in right format (e.g. S03_350622)?")

    except IOError:
        print("Could not read file {}.".format(filename))
        return None
    except Exception as err:
        print(err.args[0])
        return None
    else:
        if writeMarkerKey:
            try:
                with open(keyfilename, mode='w', newline='') as outcon:
                    cw = csv.writer(outcon)
                    cw.writerow(["TASSEL-GBSv2 marker name", "TagDigger marker name"])
                    for mk in markerkey:
                        cw.writerow(mk)
            except IOError:
                print("Could not write file {}.".format(keyfilename))
                return None
        return [namelist, seqlist]

def readTags_pyRAD(filename, toKeep = None, binaryOnly = False):
    '''Read the .alleles output from pyRAD and import tags.'''
    namelist = [] # for storing tag names
    seqlist = []  # for storing tag sequences
    theseseq = set() # set of sequences for the first marker
    allowedchars = {'A', 'C', 'G', 'T', '-', 'N'}
    linenum = 0
    
    def seqformarker(seq, m):
        '''Process a group of sequences for one marker.'''
        # trim sequences to the length of the shortest one
        seqlen = min([len(s) for s in seq])
        seq = [s[:seqlen] for s in seq]
        while any([s[-1] == '-' for s in seq]):
            seq = [s[:-1] for s in seq]
            seqlen -= 1
        # remove sequences with N's
        seq = [s for s in seq if 'N' not in s] 
        seq = sorted(set(seq)) # sort alphabetically
        nseq = len(seq)
        
        nl = [] # to add to namelist
        sl = [] # to add to seqlist

        if (nseq != 0 and not binaryOnly) or nseq == 2:
            # remove gaps from sequences
            sl = [s.replace('-','') for s in seq]
            # find variable sites for generating tag names
            alleles = [[s[i] for s in seq] for i in range(seqlen) if len(set([s[i] for s in seq])) > 1]
            alstr = [''.join([a[i] for a in alleles]) for i in range(nseq)]
            # make tag names
            nl = ['{}_{}_{}'.format(m, alstr[i], i) for i in range(nseq)]
        return [nl, sl]
    
    try:  # Read file.
        with open(filename, mode = 'r') as mycon:
            for line in mycon:
                if line[0] == '>':    # line with sequence data
                    thisseq = line.split()[1]
                    if not set(thisseq) <= allowedchars:
                        raise Exception("Character other than ACGTN- detected in sequence.")
                    theseseq.add(thisseq) # add the sequence from this line to the set for this marker
                elif line[0] == '/':  # line with marker number
                    mrkrnum = line.split()[-1][1:-1] # number for this marker
                    mrkrnum = mrkrnum.replace("|", "")
                    mrkrnum = mrkrnum.replace("*", "")
                    mrkrnum = mrkrnum.replace("-", "")
                    if toKeep == None or mrkrnum in toKeep:
                        x = seqformarker(theseseq, mrkrnum)
                        namelist.extend(x[0])  # add tag names to master list
                        seqlist.extend(x[1])   # add sequences to master list
                    theseseq = set()       # reset sequences for next marker
                else:
                    raise Exception("File not in pyRAD format.")
                linenum += 1 # increment count of which line we are on
        result = [namelist, seqlist]
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    except Exception as err:
        print("Line {}:".format(linenum))
        print(err.args[0])
        result = None
    return result
        
def readMarkerNames(filename):
    '''Read in a simple list of marker names, and use for selecting markers
       to keep from a larger list of tags.'''
    try:
        with open(filename, mode='r') as mycon:
            mylines = mycon.readlines()
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    else:
        # clear out commas, whitespace, and empty lines
        result = [x.replace(",", "").strip() for x in mylines if \
                  x.replace(",", "").strip() != ""]
    return result

def readTags_interactive():
    '''Interactively read in sequence tags.  To be used by tagdigger and
      tag manager.'''
    # read in optional file of marker names
    toKeep = None
    print('''
Do you wish to supply a list of marker names?  If provided, this list
will be used to subset the list of markers in the tag file.''')
    thischoice = ""
    while thischoice.upper() not in {'Y', 'N'}:
        thischoice = input("Y/N: ").strip()
    print("")
    if thischoice.upper() == 'Y':
        while toKeep == None:
            toKeep = readMarkerNames(input("File name: ").strip())
        # summarize file of marker names
        print('''
File contains {} marker names.'''.format(len(toKeep)))
        for i in range(min(10, len(toKeep))):
            print(toKeep[i])
        if len(toKeep) > 10:
            print('...')

    # list tag file formats:
    print('''
Available tag file formats are:
  1: UNEAK FASTA
  2: Merged tags
  3: Tags in columns
  4: Tags in rows
  5: Stacks catalog
  6: SAM file for TASSEL-GBSv2 pipeline
  7: pyRAD .alleles output
''')
    tagfunctions = {'1': readTags_UNEAK_FASTA,
                    '2': readTags_Merged,
                    '3': readTags_Columns,
                    '4': readTags_Rows,
                    '5': readTags_Stacks,
                    '6': readTags_TASSELSAM,
                    '7': readTags_pyRAD }

    # choose format and read tag file
    tags = None
    while tags == None:
        thischoice = '0'
        while thischoice not in {'1', '2', '3', '4', '5', '6', '7'}:
            thischoice = input("Enter the number of the format of your tag file: ").strip()
        if thischoice == '5':
            tagsfile = input("Enter the name of the *.catalog.tags.tsv file: ").strip()
            snpsfile = input("Enter the name of the *.catalog.snps.tsv file: ").strip()
            allelesfile = input("Enter the name of the *.catalog.alleles.tsv file: ").strip()
            binchoice = ""
            while binchoice not in {'Y', 'N'}:
                binchoice = input("Only retain binary markers? y/n: ").strip().upper()
            tags = readTags_Stacks(tagsfile, snpsfile, allelesfile, toKeep = toKeep,
                                   binaryOnly = binchoice == 'Y')
        elif thischoice == '6':
            tagfile = input("Enter the file name: ").strip()
            binchoice = ""
            while binchoice not in {'Y', 'N'}:
                binchoice = input("Only retain binary markers? y/n: ").strip().upper()
            keychoice = ""
            while keychoice not in {'Y', 'N'}:
                keychoice = input("Output a key file matching TASSEL-GBSv2 SNP names to TagDigger marker names? y/n: ").strip().upper()
            kfn = None
            if keychoice == 'Y':
                kfn = input("File name for CSV file with key: ").strip()
            print("Reading {}...".format(tagfile))
            tags = readTags_TASSELSAM(tagfile, toKeep = toKeep, binaryOnly = binchoice == 'Y',
                                      writeMarkerKey = keychoice == 'Y', keyfilename = kfn)
        elif thischoice == '7':
            tagfile = input("Enter the file name: ").strip()
            binchoice = ""
            while binchoice not in {'Y', 'N'}:
                binchoice = input("Only retain binary markers? y/n: ").strip().upper()
            print("Reading {}...".format(tagfile))
            tags = readTags_pyRAD(tagfile, toKeep = toKeep, binaryOnly = binchoice == 'Y')
        else:
            tagfile = input("Enter the file name: ").strip()
            tags = tagfunctions[thischoice](tagfile, toKeep = toKeep)
        print('')
    # summarize results
    print("{} tag sequences read.\n".format(len(tags[1])))

    return tags

def sanitizeTags(taglist):
    '''Eliminate tags that will cause problems from the list.  'taglist' is
       in the format of output from any of the readTags functions: a list of
       two elements, with the first being a list of tag names and the
       second being a list of tag sequences.'''
    assert len(taglist) == 2, "'taglist' should have two elements."
    assert len(taglist[0]) == len(taglist[1]), \
      "List of tag names should be the same as list of tag sequences."
    # Don't do a check here that tag sequences are ACGT; this is done in several
    # other places.

    print("\nSanitizing tags...")
    # Eliminate tags that are shorter versions of other tags
    sortedtags = sorted(taglist[1])
    for i in range(len(sortedtags) - 1):
        if sortedtags[i+1].startswith(sortedtags[i]) and \
           sortedtags[i] in taglist[1]: # if it is not removed yet
            myind = taglist[1].index(sortedtags[i]) # position of this tag in the list.
            tagname = taglist[0][myind]
            markername = tagname[:tagname.find("_")]
            print("Removing " + markername + " for overlap with another marker.")
            # indexes for this marker
            allind = [i for i in range(len(taglist[1])) if taglist[0][i].startswith(markername)]
            allind.sort()
            allind.reverse()
            for i in allind: # remove elements at this index
                print(taglist[0].pop(i))
                print(taglist[1].pop(i))
    return taglist


def combineReadCounts(countsdict, bckeys):
    '''Combine read counts across multiple libraries, merging samples with
       identical names.'''
    # get fastq file names
    fqfiles = sorted(bckeys.keys())
    # get total number of samples
    allsam = set()
    for f in fqfiles:
        allsam.update(bckeys[f][1])
    numsam = len(allsam)
    # get number of tags
    numtag = len(countsdict[fqfiles[0]][0])
    # set up matrix for total counts
    totcounts = [[0 for x in range(numtag)] for x in range(numsam)]
    # ordered sample names to go with the matrix
    samout = ["" for x in range(numsam)]

    # sample number index
    si = 0
    # go through read counts
    for f in fqfiles:
        for s in range(len(bckeys[f][1])):
            # if sample is already in the output
            if bckeys[f][1][s] in samout:
                # index of this sample in the matrix
                thisindex = samout.index(bckeys[f][1][s])
                # add the counts
                totcounts[thisindex] = [countsdict[f][s][i] + totcounts[thisindex][i] for \
                                        i in range(numtag)]
            # if sample is new
            else:
                # add sample name to list
                samout[si] = bckeys[f][1][s]
                # add counts to matrix
                totcounts[si] = countsdict[f][s]
                # increment new sample index
                si += 1
    return [samout, totcounts]

def writeCounts(filename, counts, samnames, tagnames):
    '''Write a CSV file of tag counts.'''
    assert len(samnames) == len(counts), "Length of samnames should be the same as length of counts."
    assert len(tagnames) == len(counts[0]), "Length of tagnames should be length of second dimension of counts."
    with open(filename, mode = 'w', newline = '') as csvfile:
        cw = csv.writer(csvfile)
        # header row with tag names
        cw.writerow([""] + tagnames)
        # rows with sample name and count data
        for s in range(len(samnames)):
            cw.writerow([samnames[s]] + counts[s])
    return

def extractMarkers(tagnames):
    '''Get marker names, allele names, and indexes from tag names.'''
    if len(tagnames) != len(set(tagnames)):
        raise Exception("Non-unique tag names found.")
    markernames = []
    alleleindex = []
    numMarkers = 0
    markernamesSort = [] # sorted list of marker names (for binary search)
    markernamesIndexSort = [] # in same order as markernamesSort, giving index in markernames
    tagindex = 0
    for t in tagnames:
        thismarker = t[:t.find('_')]
        # index of this marker in the sorted list
        thisIndSort = bisect.bisect_left(markernamesSort, thismarker)
        # add marker to list if it is new
        if thisIndSort == numMarkers or markernamesSort[thisIndSort] != thismarker:
            markernames.append(thismarker) # add to non-sorted list
            alleleindex.append([[],[]])
            # item 0 will be allele names, item 1 will be tag indices
            markernamesSort.insert(thisIndSort, thismarker) # add to sorted list
            markernamesIndexSort.insert(thisIndSort, numMarkers) # add index to sorted list
            mi = numMarkers # index for unsorted list (at end of list)
            numMarkers += 1
        else: # if marker is already in the list, get its index
            mi = markernamesIndexSort[thisIndSort]
        thisallele = t[t.rfind('_') + 1:]
        alleleindex[mi][0].append(thisallele)
        alleleindex[mi][1].append(tagindex)
        tagindex += 1
    return [markernames, alleleindex]

def writeDiploidGeno(filename, counts, samnames, tagnames):
    '''Write numeric diploid genotypes to a CSV file.'''
    assert len(samnames) == len(counts), "Length of samnames should be the same as length of counts."
    assert len(tagnames) == len(counts[0]), "Length of tagnames should be length of second dimension of counts."
    # get marker names and tag indices
    mrkr = extractMarkers(tagnames)
    nm = len(mrkr[0]) # number of markers
    ns = len(samnames) # number of samples
    try:
        if not all(set(a[0]) <= {'0', '1'} for a in mrkr[1]):
            raise Exception("All allele names must be '0' or '1'.")
        # set up genotypes
        genotypes = [['' for i in range(nm)] for j in range(ns)]
        # fill in genotypes
        for s in range(ns):
            for m in range(nm):
                counts0 = counts[s][mrkr[1][m][1][mrkr[1][m][0].index('0')]]
                counts1 = counts[s][mrkr[1][m][1][mrkr[1][m][0].index('1')]]
                if counts0 > 0 and counts1 == 0:
                    genotypes[s][m] = '0'
                if counts0 > 0 and counts1 > 0:
                    genotypes[s][m] = '1'
                if counts0 == 0 and counts1 > 0:
                    genotypes[s][m] = '2'

        with open(filename, mode = 'w', newline = '') as csvfile:
            cw = csv.writer(csvfile)
            # header row with marker names
            cw.writerow([""] + mrkr[0])
            # rows with sample name and count data
            for s in range(ns):
                cw.writerow([samnames[s]] + genotypes[s])
    except IOError:
        print("Could not write file {}.".format(filename))
    except Exception as err:
        print(err.args[0])
    return None

def set_directory_interactive():
    '''Allow user to change working directory.'''
    currentdir = os.getcwd()
    print("\nCurrent directory is:")
    print(currentdir)
    thischoice = ""
    while thischoice.upper() not in {'Y', 'N'}:
        thischoice = input("Use different directory for reading and writing files? (y/n) ").strip()
    if thischoice.upper() == 'Y':
        dirchoice = ""
        while not os.path.isdir(dirchoice):
            dirchoice = input("New directory: ")
        os.chdir(dirchoice)

    print("\nContents of current directory:")
    thisdircontents = os.listdir('.')
    for i in thisdircontents:
        print(i)
    return None

## Additional functions for barcode splitter ##
def reverseComplement(sequence):
    '''Make the reverse complement of a nucleotide sequence.'''
    x = {ord('A'): 'T', ord('C'): 'G', ord('G'): 'C', ord('T'): 'A'}
    return(sequence.translate(x)[::-1])

def build_adapter_tree(adapter, barcodes):
    '''Build a set of indexing trees, one for each barcode, for rapidly
       searching for adapter sequence at the end of a sequence.  Return
       each tree as well as a list of indices for slicing the sequences.
       Note that the sequences are reversed, in order to search from the
       end of a sequence.'''
    # common cutter adapter
    rl0 = adapter[0][0].find('^') # length of remains of restr. site
    a0 = adapter[0][0][:rl0] + adapter[0][1] # sequence to search for
    a0rev = a0[::-1] # reverse sequence
    # all potential portions of adapter sequence that could be found
    a0slices = [a0rev[i:] for i in range(len(a0rev) - rl0)]
    # index for slicing out adapter sequence of that length
    a0ind = [0 - len(a) + rl0 for a in a0slices]

    result = []
    # loop through rare cutter adapter with all barcodes
    for bi in range(len(barcodes)):
        rl1 = adapter[1][0].find('^')
        a1 = adapter[1][0][:rl1] + \
             adapter[1][1].replace('[barcode]', reverseComplement(barcodes[bi]))
        a1rev = a1[::-1]
        a1slices = [a1rev[i:] for i in range(len(a1rev) - rl1)]
        a1ind = [0 - len(a) + rl1 for a in a1slices]
        # build tree and add slicing indices
        try:
            result.append([build_sequence_tree(a0slices + a1slices, \
                                               len(a0slices + a1slices)),
                           a0ind + a1ind])
        except AssertionError:
            print("Some overlap of adapter sequence for barcode {}.".format(barcodes[bi]))
            allslices = a0slices + a1slices
            allslices.sort()
            removeIndex = []
            for ai in range(len(allslices)-1):
                if allslices[ai+1].startswith(allslices[ai]):
                    removeIndex.append(ai + 1)
                    print("Won't search for {0} at end of sequence since {1} is already being searched for.".format(allslices[ai+1][::-1], allslices[ai][::-1]))
            allslices = [allslices[ai] for ai in range(len(allslices)) if ai not in removeIndex]
            allind = [0 - len(a) + rl1 for a in allslices]
            result.append([build_sequence_tree(allslices, len(allslices)), allind])
    return result

def findAdapterSeq(sequence, adaptertree, fullsite0, fullsite1, searchstart):
    '''Find restriction site and/or adapters in a sequence, and give an
       index for slicing the string to remove them (or 999 to indicate not
       found).  Function intended to be used by barcode splitting function.
       'adaptertree' is a list of two elements, where the first is an indexing
       tree and the second is a list of indexes for each sequence to be used
       for slicing out the adapter; it is one element of the output of
       build_adapter_tree, for the appropriate barcode.
       'fullsite0' and 'fullsite1' are the full restriction sites for the common
       and rare cutters, respectively.'''
    # add functionality later for enzymes with multiple cut sites like ApeKI.

    # first see if full restriction site is present; cut short for possible
    # chimeric sequence. Returns a positive index.
    rs0 = sequence.find(fullsite0, searchstart)
    rs1 = sequence.find(fullsite1, searchstart)

    # If full restriction site note present, see if adapter is in sequence.
    # Return a negative index if so.
    if rs0 == -1 and rs1 == -1:
        adLookup = sequence_index_lookup(sequence[::-1], adaptertree[0])
        if adLookup == -1: # adapter not found
            return 999
        else:
            return adaptertree[1][adLookup]
    elif rs1 == -1: # if only the common cutter site is present
        return rs0 + len(fullsite0)
    elif rs0 == -1: # if only the rare cutter site is present
        return rs1 + len(fullsite1)
    elif rs0 < rs1: # if both are present but common cutter comes first
        return rs0 + len(fullsite0)
    else: # if both are present but rare cutter comes first
        return rs1 + len(fullsite1)


def barcodeSplitter(inputFile, barcodes, outputFiles, cutsite = 'TGCAG',
                    adapter = adapters["PstI-MspI-Hall"],
#                    maxreads = 100000):
                    maxreads=500000000):
    '''Function to split one FASTQ file into multiple files by barcode,
       removing barcode and adapter sequence.'''
    # for now only allow single cut sites
    assert set(cutsite) <= set('ACGT'), "Only ACGT cut sites allowed."
    assert all([set(bc) <= set('ACGT') for bc in barcodes]), "Found non-ACGT barcodes."
    # check adapter info
    assert len(adapter) == 2
    assert all([set(a[0]) <= set('ACGT^') for a in adapter])
    assert set(adapter[0][1]) <= set('ACGT')
    assert set(adapter[1][1]) <= set('[barcode]ACGT')

    print("Building indices for rapid searching...")

    # barcodes setup
    barlen = [len(bc) for bc in barcodes]
    barcut = combine_barcode_and_cutsite(barcodes, cutsite)
    barcuttree = build_sequence_tree(barcut, len(barcut))
    cutlen = len(cutsite)
    # adapter tree setup
    adaptertrees = build_adapter_tree(adapter, barcodes)
    # full cutsite setup
    fullsite0 = adapter[0][0].replace('^', '')
    fullsite1 = adapter[1][0].replace('^', '')

    print("Done with indexing setup.")
    print(inputFile)

    # start search through sequence
    if inputFile[-2:].lower() == 'gz':
        fqcon = gzip.open(inputFile, 'rt')
    else:
        fqcon = open(inputFile, 'r')
    readscount = 0
    barcutcount = 0
    clippedcount = 0
    lineindex = 0
    outcons = [open(outfile, mode = 'w') for outfile in outputFiles]
    try:
        for line in fqcon:
            if lineindex % 4 == 0:
                comment1 = line.strip()
            if lineindex % 4 == 1:
                sequence = line.strip().upper()
            if lineindex % 4 == 2:
                comment2 = line.strip()
            if lineindex % 4 == 3:
                readscount += 1
                quality = line.strip()
                barindex = sequence_index_lookup(sequence, barcuttree)
                if barindex > -1: # if it matches a barcode
                    barcutcount += 1
                    # get indexes for clipping sequence
                    slice1 = barlen[barindex]
                    slice2 = findAdapterSeq(sequence, adaptertrees[barindex],
                                            fullsite0, fullsite1, slice1+cutlen)
                    if slice2 == 999: # if no trimming from right side is necessary
                        slice2 = len(sequence)
                    else:
                        clippedcount += 1
                    # write clipped output
                    outcons[barindex].write(comment1 + barcodes[barindex] + '\n')
                    outcons[barindex].write(sequence[slice1:slice2] + '\n')
                    if comment2 == '+':
                        outcons[barindex].write('+\n')
                    else:
                        outcons[barindex].write(comment1 + barcodes[barindex] + '\n')
                    outcons[barindex].write(quality[slice1:slice2] + '\n')
                if readscount % 1000000 == 0:
                    print(inputFile)
                if readscount % 50000 == 0:
                    print("Reads: {0} With barcode and cut site: {1} Clipped on 3' end: {2}".format(readscount, barcutcount, clippedcount))
                if readscount >= maxreads:
                    break
            lineindex += 1
    finally:
        fqcon.close()
        for o in outcons:
            o.close()
    return None

def writeMD5sums(filelist, outfile):
    '''Write a CSV of file names and MD5 checksums, given a list of file names.'''
    maxfilelen = max([len(f) for f in filelist])
    with open(outfile, mode='w', newline = '') as csvcon:
        cw = csv.writer(csvcon)
        cw.writerow(["File name", "MD5 sum"])
        for f in filelist:
            m = hashlib.md5() # variable to contain MD5 sum
            with open(f, 'rb') as fqcon:
                while True:
                    chunk = fqcon.read(50 * 1048576) # read 50 Mb at a time
                    if chunk == b'': # end of file
                        break
                    m.update(chunk)
            cw.writerow([f, m.hexdigest()])
            print("{:>{width}} {}".format(f, m.hexdigest(), width=maxfilelen))
    return None

## Additional functions for tag manager ##
def exportFasta(filename, namelist, seqlist):
    '''Function to take a list of tags and export to a FASTA file for use in
       an alignment program such as Bowtie2 or BLAST.  For multiple tags
       belonging to one marker, include ambiguous nucleotide codes.'''
    assert len(namelist) == len(seqlist), "List of marker names and list of tag sequences should be same length."
    assert all([set(t) <= set('ACGT') for t in seqlist]), "Tag sequences need to be ACGT."

    # get indexes of alleles by marker
    markerindex = extractMarkers(namelist)

    try:
        with open(filename, mode='w') as mycon:
            for mi in range(len(markerindex[0])): # loop through markers
                markername = markerindex[0][mi]
                # forbid spaces in marker names (will screw with Bowtie2)
                if ' ' in markername:
                    raise Exception("{}: Marker names cannot contain spaces.".format(markername))
                # write comment line with marker name
                mycon.write('>' + markername + '\n')
                # get list of tags for this marker
                mtags = [seqlist[i] for i in markerindex[1][mi][1]]
                if len(mtags) == 1: # if non-variable
                    mycon.write(mtags[0] + '\n')
                else:
                    # find variable sites
                    ctags = compareTags(mtags, trim = False)
                    # write the first non-variable portion of the tag
                    mycon.write(mtags[0][:ctags[0][0]])
                    # cycle through variable sites and the following sequence
                    for c in range(len(ctags)):
                        thisset = set(ctags[c][1])
                        thisset.discard('N')
                        mycon.write(IUPAC_codes[frozenset(thisset)])
                        if c == len(ctags)-1:
                            mycon.write(mtags[0][ctags[c][0] + 1:])
                        else:
                            mycon.write(mtags[0][ctags[c][0] + 1:ctags[c+1][0]])
                    mycon.write('\n')

    except IOError:
        print("Could not write file {}.".format(filename))
    except Exception as err:
        print(err.args[0])
        os.remove(filename)
    return None

def varSitesByMarker(namelist, seqlist):
    '''Across a set of tags as imported by a readTags function, make a dictionary of
       the positions of all variable sites for each marker.'''
    markers = extractMarkers(namelist)
    nMarkers = len(markers[0]) # number of markers
    outDict = dict() # dictionary output
    # loop through markers
    for mi in range(nMarkers):
        tagindex = markers[1][mi][1] # index of tags for this marker in the list of sequences
        # comparison of tags
        comp = compareTags([seqlist[i] for i in tagindex])
        # add variable sites to dictionary
        outDict[markers[0][mi]] = [j[0] for j in comp]
    return outDict

def readSAM(filename, varDict = None):
    '''Read in a SAM file of alignment information for markers.  Return a
       dictionary, with marker names as the keys, and tuples of reference
       sequences, positions, and quality scores as the values.
       varDict is output from varSitesByMarker.  If not None, the last item
       in the tuple for each marker will be a list of variable sites.'''
    result = dict()
    try:
        with open(filename, mode='r') as mycon:
            for line in mycon:
                if line[0] == '@': # skip headers
                    continue
                mycolumns = line.split()
                myflags = int(mycolumns[1])
                # skip if no alignment (4 flag)
                if myflags - 4 in {0, 1, 2, 8, 16, 32, 64, 128}:
                    continue
                # calculate actual variable positions if desired
                if varDict != None:
                    if myflags - 16 in {0, 1, 2, 8, 32, 64, 128}: # bottom strand
                        taglen = len(mycolumns[9]) # length of sequence
                        diff = [int(mycolumns[3]) + taglen - 1 - i for i in varDict[mycolumns[0]]]
                    else: # top strand
                        diff = [int(mycolumns[3]) + i for i in varDict[mycolumns[0]]]
                else:
                    diff = None

                # column 1 is marker name, 3 is chromosome, 4 is position, 5 is quality
                if varDict != None:
                    result[mycolumns[0]] = (mycolumns[2], mycolumns[3], mycolumns[4], diff)
                else:
                    result[mycolumns[0]] = (mycolumns[2], mycolumns[3], mycolumns[4])
        return result
    except IOError:
        print("Could not read file {}.".format(filename))
        return None
    except IndexError:
        print("File {} in wrong format.".format(filename))
        return None

def mergeTags(tags):
    '''Given a list of sequences, produce an output string in "merged" format,
       with square brackets surrounding the variable region and a forward
       slash separating the two variants.'''
    taglen = [len(t) for t in tags]
    ntags = len(tags)
    longesttag = [tags[i] for i in range(ntags) if taglen[i] == max(taglen)][0]
    x = compareTags(tags, trim = False) # has assert to check ACGT
    assert len(x) > 0, "All tags in set are identical."
    variantpositions = [y[0] for y in x]
    minvar = min(variantpositions)
    maxvar = max(variantpositions) if len(set(taglen)) == 1 else max(taglen) - 1
    # get invariant portions of tags
    invarstart = longesttag[:minvar]
    invarend = longesttag[maxvar + 1:]
    # get variant portions
    var = [t[minvar:(maxvar+1)] if len(t) > maxvar else t[minvar:] for t in tags]
    return invarstart + '[' + '/'.join(var) + ']' + invarend

def mergedTagList(tags):
    '''Given a list of tag names and tags such as those output by the
       readTags functions, output a list of marker names and a list of
       merged tag strings (generated by mergeTags).'''
    markers = extractMarkers(tags[0]) # marker names, tag indices
    nM = len(markers[0]) # number of markers
    try:
        # check that each marker has two alleles
        if not all([len(m[1]) > 1 for m in markers[1]]):
            raise Exception("Each marker needs multiple tags.")
        mergedStrings = ["" for i in range(nM)] # list to contain merged strings
        for i in range(nM):
            # get indices of tags for this marker, sorted by allele name
            tagIndices = [ti for (al, ti) in sorted(zip(markers[1][i][0], markers[1][i][1]))]
            mergedStrings[i] = mergeTags([tags[1][ti] for ti in tagIndices])
        return [markers[0], mergedStrings]
    except Exception as err:
        print(err.args[0])
        return None

def exportFasta2(filename, markernames, mergedstrings):
    '''Take the merged tag strings produced by mergedTagList and write a FASTA
       file for alignment with software such as Bowtie2 or BLAST.  Use IUPAC
       codes for variable nucleotides.
       Deprecated function that only allows for two tags per marker'''
    assert len(markernames) == len(mergedstrings), \
           "Must have same number of marker names and merged strings."
    try:
        with open(filename, mode='w') as mycon:
            for i in range(len(markernames)):
                if ' ' in markernames[i]:
                    raise Exception("{}: Marker names cannot contain spaces.".format(markernames[i]))
                thisstring = mergedstrings[i]
                if not set(thisstring) <= set('[/]ACGT'):
                    raise Exception("{}: Unexpected character in merged string.".format(markernames[i]))
                if not set('[/]') < set(thisstring):
                    raise Exception("{}: Square brackets and slash not found.".format(markernames[i]))
                mycon.write('>' + markernames[i] + '\n') # comment line
                p1 = thisstring.find('[')
                p2 = thisstring.find('/')
                p3 = thisstring.find(']')
                if p1 > p2 or p2 > p3:
                    raise Exception("{}: Square brackets and slash in wrong order.".format(markernames[i]))
                mycon.write(thisstring[:p1]) # first non-variable section
                var1 = thisstring[p1+1:p2]
                var2 = thisstring[p2+1:p3]
                if len(var1) != len(var2):
                    raise Exception("{}: Variable regions are of different lengths.".format(markernames[i]))
                for j in range(len(var1)): # variable section
                    mycon.write(IUPAC_codes[frozenset({var1[j], var2[j]})])
                mycon.write(thisstring[p3+1:] + '\n') # last non-variable section
    except IOError:
        print("Could not write file {}.".format(filename))
    except Exception as err:
        print(err.args[0])
        os.remove(filename)
    return None

def readTabularData(filename, markerDict = None, ignoreSeq = False):
    '''Read in a table of extra columns to include in tag database.
       Optionally, markerDict is a dictionary for converting marker
       names.  You can turn two lists into dictionary for markerDict
       using dict(zip(list1, list2)).'''
    try:
        with open(filename, 'r', newline='') as mycon:
            mycsv = csv.reader(mycon)
            rowcount = 0
            dataDict = dict() # marker names are keys, items are rows
            for row in mycsv:
                if rowcount == 0:
                    if "Marker name" not in row:
                        raise Exception("Need a 'Marker name' column header.")
                    mi = row.index("Marker name")
                    headers = row
                    headers.pop(mi)
                    if ignoreSeq: # if we want to ignore any column header containing sequence
                        si = row.index("Tag sequence")
                        headers.pop(si)
                else:
                    thismarker = row.pop(mi)
                    if markerDict != None and thismarker in markerDict.keys():
                        thismarker = markerDict[thismarker]
                    if ignoreSeq:
                        row.pop(si)
                    if thismarker in dataDict:
                        # consolidate two rows into one
                        dataDict[thismarker] = [row[i] if row[i].strip() != "" else \
                                                dataDict[thismarker][i] for i in range(len(row))]
                    else:
                        dataDict[thismarker] = row
                rowcount += 1
        return [headers, dataDict]
    except IOError:
        print("Could not read file {}.".format(filename))
        return None
    except Exception as err:
        print(err.args[0])
        return None

def writeMarkerDatabase(filename, markernames, mergedseq, extracollist):
    '''Write a CSV giving sequence and other information for markers.  First
       column is marker names.  Second column is sequence in merged format.
       Subsequent columns are described by extracollist.  Each item in
       extracollist has two items, the first being a list of column headers,
       and the second being a dictionary with marker names as the key, and the
       item being a list of data, in order, for the extra columns.'''
    assert isinstance(extracollist, list), "extracollist must be a list (empty if not needed)."
    assert all([len(l) == 2 for l in extracollist]), "Each item in extracollist needs two components."
    assert all([isinstance(l[1], dict) for l in extracollist]), "extracollist needs dictionaries."
    try:
        with open(filename, 'w', newline = '') as mycon:
            cw = csv.writer(mycon)
            headernames = ['Marker name', 'Tag sequence']
            nExtraCol = [] # number of extra columns, for each dictionary in extracollist
            nList = len(extracollist)
            for l in extracollist: # add header names for extra columns
                headernames.extend(l[0])
                nExtraCol.append(len(l[0]))
            cw.writerow(headernames)
            for i in range(len(markernames)): # loop through markers and write rows
                m = markernames[i] # this marker
                thisrow = [m, mergedseq[i]]
                for j in range(nList):
                    d = extracollist[j][1] # the dictionary
                    if m not in d.keys():
                        thisrow.extend(["" for i in range(nExtraCol[j])])
                    else:
                        thisrow.extend(d[m])
                cw.writerow(thisrow)
    except IOError:
        print("Could not write file {}.".format(filename))
    return None

def readMarkerDatabase(filename):
    '''Read in a CSV file containing a SNP marker database as produced by
       writeMarkerDatabase.  Read in marker names, tag sequence, and any
       extra columns.'''
    print("Reading data...")
    try:
        tags = readTags_Merged(filename, allowDuplicates=True) # read tag sequence
        if tags == None:
            raise IOError # should happen if file is not readable
        addData = readTabularData(filename, ignoreSeq = True) # read other columns
        if addData == None:
            raise IOError
        return [tags, addData]
    except IOError:
#        print("Could not read file {}.".format(filename))
        return None
    except Exception as err:
        print(err.args[0])
        return None

def lookupMarkerByTag(tagNamesSort, tagSeqSort, queryTags, allowDiffLengths = False):
    '''For a set of tags corresponding to one marker for the query, retrieve a
       set of markers with matching tags.  tagSeqSort is sorted list of tag sequences,
       tagNamesSort is a list of tag names in the same order (markername_allelename),
       and queryTags is a list of tag sequences for one marker.
       If allowDiffLengths = True, the tags don't need to be the same length to match, 
       but instead one can start with the other.  The tag in queryTags will be ignored
       if more than one tag in tagSeqSort starts with it.  The function will however
       accomodate multiple identical tags in tagSeqSort.'''
    markersOut = set()
    nTags = len(tagSeqSort)
    assert nTags == len(tagNamesSort), "tagNamesSort and tagSeqSort not same length"
    for tag in queryTags:
        markerfound = False
        sortind = bisect.bisect_left(tagSeqSort, tag)
        if sortind < nTags and tag == tagSeqSort[sortind]:
            thistagname = tagNamesSort[sortind]
            markerfound = True
        elif allowDiffLengths:
            if sortind > 0 and tag.startswith(tagSeqSort[sortind - 1]):
                thistagname = tagNamesSort[sortind - 1]
                markerfound = True
                # adjust the matching index to check for multiples later
                sortind -= 1
                while sortind > 0 and tagSeqSort[sortind] == tagSeqSort[sortind - 1]:
                    sortind -= 1
            if sortind < nTags and tagSeqSort[sortind].startswith(tag) and not \
              (sortind < nTags - 1 and tagSeqSort[sortind] != tagSeqSort[sortind + 1] and tagSeqSort[sortind + 1].startswith(tag)):
                thistagname = tagNamesSort[sortind]
                markerfound = True
        if markerfound:
            thismarkername = thistagname[:thistagname.find('_')]
            markersOut.add(thismarkername)
            # add more if it matches multiples
            while sortind < len(tagSeqSort) - 1 and (tagSeqSort[sortind] == tagSeqSort[sortind + 1] or \
              (allowDiffLengths and tagSeqSort[sortind + 1].startswith(tag))):
                thistagname = tagNamesSort[sortind + 1]
                thismarkername = thistagname[:thistagname.find('_')]
                markersOut.add(thismarkername)
                sortind += 1
            while allowDiffLengths and sortind > 0 and tag.startswith(tagSeqSort[sortind - 1]):
                thistagname = tagNamesSort[sortind - 1]
                thismarkername = thistagname[:thistagname.find('_')]
                markersOut.add(thismarkername)
                sortind -= 1
                
    return markersOut
    
def sortTagsBySeq(tags):
    '''Take tags output by one of the readTags functions, and return tags and tag
    names in the same format but sorted by sequence (for binary search).'''
    seqs, names = zip(*sorted(zip(tags[1], tags[0])))
    return [names, seqs]

def compareTagSets(oldtags, newtags, perfectMatch = False, allowDiffLengths = True):
    '''Compare two sets of tags, in the format output by the readTags function.
       Return a dictionary where the keys include all marker names from newtags,
       and the items are the marker names from oldtags (if found).
       perfectMatch = does each tag in the new set need to have a match in the old set and vice versa?
       allowDiffLengths = can a new tag be a shorter version of an old tag and vice versa?'''
    oldtags_sort = sortTagsBySeq(oldtags) # old tag sequences sorted for binary searching
    oldmarkers = extractMarkers(oldtags_sort[0]) # get marker indices in tag lists
    newmarkers = extractMarkers(newtags[0])
    NnewMarkers = len(newmarkers[0]) # number of new markers
    Nnewtags = len(newtags[0]) # number of new tags (normally twice the number of markers)
    resultDict = dict.fromkeys(set(newmarkers[0])) # output dictionary

    # build structures for binary searching
    oldmarkerIndSort = [x[1] for x in sorted(zip(oldmarkers[0], range(len(oldmarkers[0]))))]
    oldmarkerSort = sorted(oldmarkers[0])

    for m in range(NnewMarkers): # loop through markers
        thismarker = newmarkers[0][m] # marker name
        resultDict[thismarker] = []
        # tag sequences for this marker
        theseseq = [newtags[1][i] for i in newmarkers[1][m][1]]
        mrkrmatch = lookupMarkerByTag(oldtags_sort[0], oldtags_sort[1], theseseq, allowDiffLengths = allowDiffLengths)
        if perfectMatch and len(mrkrmatch) == 1:
            oldmarker = mrkrmatch.pop()
            oi = oldmarkerIndSort[bisect.bisect_left(oldmarkerSort, oldmarker)]
            oldseq = [oldtags_sort[1][i] for i in oldmarkers[1][oi][1]]
            if allowDiffLengths:
                minlen = min([len(seq) for seq in theseseq + oldseq])
                oldseq = [seq[:minlen] for seq in oldseq]
                theseseq = [seq[:minlen] for seq in theseseq]
            if set(oldseq) == set(theseseq):
                resultDict[thismarker].append(oldmarker)
        elif not perfectMatch:
            resultDict[thismarker].extend(mrkrmatch)
    return resultDict
    
def consolidateTagSets(oldtags, newtags = None, allowDiffLengths = True,
                       prefix = "Mrkr", numdig = 7, startnumnew = 1):
    '''Allowing for multiple tags per marker, and for tag sets that do not completely overlap
       between old and new markers, make a consolidated set of tags and marker names.  For use
       by Tag Manager when there are multiple alleles.
       prefix, numdig, and startnumnew are for indicating how new markers should be named.'''
    oldtags_sort = sortTagsBySeq(oldtags) # old tag sequences sorted for binary searching
    oldmarkers = extractMarkers(oldtags[0]) # get marker indices in tag lists
    oldmarkerSort = sorted(oldmarkers[0])
    oldmarkerIndSort = [x[1] for x in sorted(zip(oldmarkers[0], range(len(oldmarkers[0]))))]
   
    # first, find any to consolidate in old tags
    oldtags_consolidated = [[], []]
    markerMatchDict = dict() # keys are markers, items are lists of markers merged into them
    dupmarkers = set() # markers that have been consolidated
    for mi in range(len(oldmarkers[0])):
        thismarker = oldmarkers[0][mi]
        if thismarker in dupmarkers:
            continue # skip if this one was already merged into another marker
        theseseq = [oldtags[1][i] for i in oldmarkers[1][mi][1]]
        thislookup = lookupMarkerByTag(oldtags_sort[0], oldtags_sort[1], theseseq, allowDiffLengths = allowDiffLengths)
        assert thismarker in thislookup, "Marker {} not found in lookup".format(thismarker)
        thislookup.remove(thismarker)
        for mrkr in thislookup: # loop through additional markers that match this one
            dupmarkers.add(mrkr)
            di = oldmarkerIndSort[bisect.bisect_left(oldmarkerSort, mrkr)]
            assert oldmarkers[0][di] == mrkr, "Duplicate marker mismatch at {} {}".format(thismarker, mrkr)
            seqtoadd = [oldtags[1][i] for i in oldmarkers[1][di][1] if oldtags[1][i] not in theseseq]
            if allowDiffLengths: # if one tag is shorter version of other, keep the longer version
                for sNew in seqtoadd:
                    for sOldI in range(len(theseseq)):
                        sOld = theseseq[sOldI]
                        if sOld.startswith(sNew): # don't need to add is new is just shorter version of old
                            if sNew in seqtoadd: # check that the new tag hasn't already been removed
                                seqtoadd.remove(sNew)
                        if sNew.startswith(sOld):
                            theseseq[sOldI] = sNew
                            seqtoadd.remove(sNew)
            theseseq.extend(seqtoadd)
        # add sequences and tag names to list
        thiscompare = compareTags(theseseq)
        nseq = len(theseseq)
        nSNPs = len(thiscompare)
        allelenames = ["".join([thiscompare[i][1][t] for i in range(nSNPs)]) for t in range(nseq)]
        tagnames = ["{}_{}_{}".format(thismarker, allelenames[i], i) for i in range(nseq)]
        oldtags_consolidated[0].extend(tagnames)
        oldtags_consolidated[1].extend(theseseq)
        markerMatchDict[thismarker] = sorted(thislookup)
    print("{} markers consolidated into {} markers".format(len(oldmarkers[0]), len(markerMatchDict)))

    if newtags == None:
        tagsOut = oldtags_consolidated
    else: 
        newtemp = consolidateTagSets(newtags, newtags = None, allowDiffLengths = allowDiffLengths)
        newtags_consolidated = newtemp[0]
        newtags_sort = sortTagsBySeq(newtags_consolidated)
        newmarkers = extractMarkers(newtags_consolidated[0])
        newmarkerSort = sorted(newmarkers[0])
        newmarkerIndSort = [x[1] for x in sorted(zip(newmarkers[0], range(len(newmarkers[0]))))]
        oldmarkers = extractMarkers(oldtags_consolidated[0])
        tagsOut = [[], []] # tags to output (same format as readTags functions)
        newmarkers_matched = set()
        for mi in range(len(oldmarkers[0])):
            thismarker = oldmarkers[0][mi]
            theseseq = [oldtags_consolidated[1][i] for i in oldmarkers[1][mi][1]]
            thislookup = lookupMarkerByTag(newtags_sort[0], newtags_sort[1], theseseq, allowDiffLengths = allowDiffLengths)
            for mrkr in thislookup:
                newmarkers_matched.add(mrkr)
                di = newmarkerIndSort[bisect.bisect_left(newmarkerSort, mrkr)]
                seqtoadd = [newtags_consolidated[1][i] for i in newmarkers[1][di][1] if newtags_consolidated[1][i] not in theseseq]
                if allowDiffLengths: # if one tag is shorter version of other, keep the longer version
                    for sNew in seqtoadd:
                        for sOldI in range(len(theseseq)):
                            sOld = theseseq[sOldI]
                            if sOld.startswith(sNew): # don't need to add is new is just shorter version of old
                                if sNew in seqtoadd: # check that the new tag hasn't already been removed
                                    seqtoadd.remove(sNew)
                            if sNew.startswith(sOld):
                                theseseq[sOldI] = sNew
                                seqtoadd.remove(sNew)
                theseseq.extend(seqtoadd)
                markerMatchDict[thismarker].append(mrkr)
                markerMatchDict[thismarker].extend(newtemp[1][mrkr]) # new markers consolidated into this new marker
            # add sequences and tag names to list
            thiscompare = compareTags(theseseq)
            nseq = len(theseseq)
            nSNPs = len(thiscompare)
            allelenames = ["".join([thiscompare[i][1][t] for i in range(nSNPs)]) for t in range(nseq)]
            tagnames = ["{}_{}_{}".format(thismarker, allelenames[i], i) for i in range(nseq)]
            tagsOut[0].extend(tagnames)
            tagsOut[1].extend(theseseq)
        # add any new markers that weren't a match for an old marker.  Give them new names.
        for mi in range(len(newmarkers[0])):
            thismarker = newmarkers[0][mi]
            if thismarker in newmarkers_matched:
                continue # skip if matched
            theseseq = [newtags_consolidated[1][i] for i in newmarkers[1][mi][1]]
            newmrkrname = "{}{:0{width}}".format(prefix, startnumnew, width = numdig)
            startnumnew += 1
            thesetagnames = [newtags_consolidated[0][i].replace(thismarker, newmrkrname) \
              for i in newmarkers[1][mi][1]] 
            tagsOut[0].extend(thesetagnames)
            tagsOut[1].extend(theseseq)
            markerMatchDict[newmrkrname] = [thismarker]
            markerMatchDict[newmrkrname].extend(newtemp[1][thismarker])
        print("{} markers consolidated into {} markers".format(len(oldmarkers[0]) + len(newmarkers[0]), \
                                                               len(markerMatchDict)))
    return [tagsOut, markerMatchDict]

def allColumns(extracollist):
    '''Combine all column names into one list.'''
    output = []
    for ecol in extracollist:
        output.extend(ecol[0])
    return output

def consolidateExtraCols(extracollist):
    '''Consolidate data from columns with the same column header.  Data from
       earlier items in extracollist is overwritten by data from later columns
       in extracollist.  See writeMarkerDatabase for structure of extracollist.'''
    ac = allColumns(extracollist) # all column header names
    while len(set(ac)) < len(ac): # test if some columns still need to be consolidated
        nList = len(extracollist) # number of tables to combine
        for j in range(0, nList - 1): # iterate through all pairs
            for k in range(j + 1, nList):
                # test if there are any overlapping columns between this pair
                if len(set(extracollist[j][0]) & set(extracollist[k][0])) > 0:
                    # set up new tables
                    newJ = [[e for e in extracollist[j][0] if e not in extracollist[k][0]], dict()]
                    newK = [[e for e in extracollist[k][0] if e not in extracollist[j][0]], dict()]
                    combined = [[e for e in extracollist[j][0] if e in extracollist[k][0]], dict()]
                    # indexes of columns from table J
                    newJindex = [extracollist[j][0].index(i) for i in newJ[0]]
                    combJindex = [extracollist[j][0].index(i) for i in combined[0]]
                    # update tables marker by marker
                    for jMrkr in extracollist[j][1].keys():
                        newJ[1][jMrkr] = [extracollist[j][1][jMrkr][i] for i in newJindex]
                        combined[1][jMrkr] = [extracollist[j][1][jMrkr][i] for i in combJindex]
                    # indexes of columns from table K
                    newKindex = [extracollist[k][0].index(i) for i in newK[0]]
                    combKindex = [extracollist[k][0].index(i) for i in combined[0]]
                    for kMrkr in extracollist[k][1].keys():
                        newK[1][kMrkr] = [extracollist[k][1][kMrkr][i] for i in newKindex]
                        combined[1][kMrkr] = [extracollist[k][1][kMrkr][i] for i in combKindex]
                    # update list of tables
                    extracollist[j] = newJ
                    extracollist[k] = newK
                    extracollist.append(combined)
        # cleanup
        extracollist = [ecol for ecol in extracollist if len(ecol[0]) > 0]
        ac = allColumns(extracollist)

    return extracollist
