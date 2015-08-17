# Functions and restriction enzyme data for the Python program TagDigger,
# created by Lindsay V. Clark.

import gzip
import csv

# dictionary of restriction enzyme cut sites, including only what will show
# up after the barcode.
enzymes = {'ApeKI': 'CWGC', 'EcoT22I': 'TGCAT', 'NcoI': 'CATGG',
           'NsiI': 'TGCAT', 'PstI': 'TGCAG', 'SbfI': 'TGCAGG'}

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
    return(tree_recursive(seqlist))

def sequence_index_lookup(sequence, seqtree):
    '''Lookup a sequence in an index tree, and get the index.  Return -1 if
       the sequence is not in the index tree.'''
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
                    maxreads=500000000):
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
            
    # Write matrix to file
#    with open(output, 'w') as outcon:
#        outcon.writelines('\t'.join(str(j) for j in i) + '\n' for i in mycounts)
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

def readBarcodeKeyfile(filename):
    '''Read in a csv file containing file names, barcodes, and sample names,
       and output a dictionary containing this information.'''
    try:
        with open(filename, 'r', newline='') as mycon:
            mycsv = csv.reader(mycon)
            rowcount = 0
            result = dict()
            for row in mycsv:
                if rowcount == 0:
                    # read header row
                    fi = row.index("File")
                    bi = row.index("Barcode")
                    si = row.index("Sample")
                else:
                    f = row[fi].strip() # file name
                    b = row[bi].strip().upper() # barcode
                    s = row[si].strip() # sample
                    # skip blank lines
                    if f == "" and b == "" and s == "":
                        next
                    # check row contents
                    if f == "":
                        raise Exception("Blank cell found where file name should be in row {}.".format(rowcount + 1))
                    if b == "":
                        raise Exception("Blank cell found where barcode should be in row {}.".format(rowcount + 1))
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
    except IOError:
        print("Could not read file {}.".format(filename))
        result = None
    except ValueError:
        print("File header needed containing 'File', 'Barcode', and 'Sample'.")
        result = None
    except Exception as err:
        print(err.args[0])
        result = None
    return result

def compareTags(taglist):
    '''Find SNP alleles within a set of sequence tags representing the
       same locus.'''
    assert type(taglist) is list, "taglist must be list."
    assert all([set(t) <= set('ATCG') for t in taglist]), \
           "taglist must be a list of ACGT strings."
    # make sure all tags are same length
    if len(set([len(t) for t in taglist])) > 1:
        minlen = min([len(t) for t in taglist])
        taglist = [tag[:minlen] for tag in taglist]
    # get index(es) and nucleotides for variable sites
    result =[(i, [t[i] for t in taglist]) for i in range(len(taglist[0])) \
                 if len(set([t[i] for t in taglist])) > 1]
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
                    seqlist.extend([seq1, seq2])
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

def readTags_Merged(filename, toKeep = None):
    '''Read a CSV of sequence tags, where each row has a marker name and a
       tag sequence with polymorphic nucleotides in square brackets separated
       by a forward slash.'''
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
                        rowcount +=1
                        continue # skip this marker if it is not in the list of ones we want.
                    # find delimiting characters within sequence
                    p1 = row[ti].find('[')
                    p2 = row[ti].find('/')
                    p3 = row[ti].find(']')
                    # extract the two tag sequences
                    tag1 = row[ti][:p1] + row[ti][p1+1:p2] + row[ti][p3+1:]
                    tag1 = tag1.upper().strip()
                    tag2 = row[ti][:p1] + row[ti][p2+1:p3] + row[ti][p3+1:]
                    tag2 = tag2.upper().strip()
                    if (tag1 in seqlist) or (tag2 in seqlist):
                        raise Exception("Non-unique sequence found: line {0}.".format(rowcount+1))
                    seqlist.append(tag1)
                    seqlist.append(tag2)
                    if not set(tag1 + tag2) <= set('ACGT'):
                        raise Exception("Tag sequence not formatted correctly in row {}.".format(rowcount+1))
                    # generate the two allele names
                    namelist.append(mname + '_' + row[ti][p1+1:p2] + '_0')
                    namelist.append(mname + '_' + row[ti][p2+1:p3] + '_1')
                rowcount += 1
        result = [namelist, seqlist]
    except IOError:
        print("File {} not readable.".format(filename))
        result = None
    except Exception as err:
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
    
