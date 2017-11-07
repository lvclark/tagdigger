#!/usr/bin/env python3

# script for running TagDigger in a non-interactive (scripted) way

import tagdigger_fun
import argparse
import os

# set up arguments to the script
parser = argparse.ArgumentParser(description="TagDigger v. 1.1 command line script by Lindsay V. Clark")
parser.add_argument('-e', '--enzyme', help='Restriction enzyme name', choices = sorted(tagdigger_fun.enzymes.keys()))
parser.add_argument('-c', '--cutsite', help='Restriction cut site sequence expected in sequencing reads')
parser.add_argument('-w', '--directory', help='Working directory')

# Tag files in various formats
parser.add_argument('--UNEAKtags', help = 'File name for tags in UNEAK format')
parser.add_argument('--MergedTags', help = 'File name for tags in merged format')
parser.add_argument('--ColumnTags', help = 'File name for tags in column format')
parser.add_argument('--RowTags', help = 'File name for tags in row format.')
parser.add_argument('--StacksTags', help = 'File name for Stacks tags.tsv file.')
parser.add_argument('--StacksSnps', help = 'File name for Stacks snps.tsv file.')
parser.add_argument('--StacksAlleles', help = 'File name for Stacks alleles.tsv file.')
parser.add_argument('--TASSELSAM', help = 'File name for TASSEL SAM file')
parser.add_argument('--pyRADalleles', help = 'File name for pyRAD .alleles file.')
parser.add_argument('-k', '--tokeep', help = 'File name listing tags to keep')
parser.add_argument('--binaryOnly', help = "'T' to retain only binary markers; 'F' to retain all markers.", default = 'F', choices = ['T', 'F'])
parser.add_argument('--TASSELkeyFile', help = 'File name to output for key to TASSEL SNP names')

parser.add_argument('-b', '--barcodefile', help = 'Name of barcode key file',
                    required = True)
parser.add_argument('-o', '--outputcounts', help = 'Output file name for read counts',
                    required = True)
parser.add_argument('-g', '--outputgen', help = 'Output file name for numeric genotypes')

args = parser.parse_args()

# check restriction cut site and set that variable
if args.enzyme == None and args.cutsite == None:
    raise Exception("Need either restriction enzyme name or cutsite sequence.  Use '-e None' if no restriction site is present in reads.")
if args.enzyme != None and args.cutsite != None:
    cutsite = args.cutsite.upper()
    if cutsite != tagdigger_fun.enzymes[args.enzyme]:
        raise Exception("Restriction enzyme name and cutsite do not match.  Note that only one of these two arguments is required.")
elif args.enzyme != None:
    cutsite = tagdigger_fun.enzymes[args.enzyme]
elif args.cutsite != None:
    cutsite = args.cutsite.upper()
    if not set(cutsite) <= set('ACGTRYSWKMBDHVN'):
        raise Exception("Cut site contains unexpected characters.")

# move to working directory
if args.directory != None:
    if not os.path.isdir(args.directory):
        raise Exception("Directory {} not found".format(args.directory))
    os.chdir(args.directory)

# determine tag file format
tagformat = [args.UNEAKtags != None, args.MergedTags != None,
             args.ColumnTags != None, args.RowTags != None,
             args.StacksTags != None, args.StacksSnps != None,
             args.StacksAlleles != None, args.TASSELSAM != None,
             args.pyRADalleles != None]
# check all three files for stacks
if tagformat[4] or tagformat[5] or tagformat[6]:
    if not (tagformat[4] and tagformat[5] and tagformat[6]):
        raise Exception("Need all three files for Stacks format.")
del tagformat[5:7] # remove redundant Stacks results
if sum(tagformat) != 1:
    raise Exception('Exactly one tag format required.')
# get list of markers to keep, if provided
if args.tokeep == None:
    toKeep = None
else:
    toKeep = tagdigger_fun.readMarkerNames(args.tokeep)
    if toKeep == None:
        raise Exception("Problem reading marker names to keep.")
# only retain binary markers?
binaryOnly = args.binaryOnly == 'T'
# read in tags
if tagformat[0]:
    tags = tagdigger_fun.readTags_UNEAK_FASTA(args.UNEAKtags, toKeep = toKeep)
if tagformat[1]:
    tags = tagdigger_fun.readTags_Merged(args.MergedTags, toKeep = toKeep)
if tagformat[2]:
    tags = tagdigger_fun.readTags_Columns(args.ColumnTags, toKeep = toKeep)
if tagformat[3]:
    tags = tagdigger_fun.readTags_Rows(args.RowTags, toKeep = toKeep)
if tagformat[4]:
    tags = tagdigger_fun.readTags_Stacks(args.StacksTags, args.StacksSnps, 
                                         args.StacksAlleles, 
                                         toKeep = toKeep, binaryOnly = binaryOnly)
if tagformat[5]:
    writeMarkerKey = args.TASSELkeyFile != None
    tags = tagdigger_fun.readTags_TASSELSAM(args.TASSELSAM, toKeep = toKeep, 
                              binaryOnly = binaryOnly, 
                              writeMarkerKey = writeMarkerKey,
                              keyfilename = args.TASSELkeyFile)
if tagformat[6]:
    tags = tagdigger_fun.readTags_pyRAD(args.pyRADalleles, toKeep = toKeep, 
                                        binaryOnly = binaryOnly)
if tags == None:
    raise Exception("Problem reading tags.")
tags = tagdigger_fun.sanitizeTags(tags)
    
# read barcode file
bckeys = tagdigger_fun.readBarcodeKeyfile(args.barcodefile)
if bckeys == None:
    raise Exception("Problem reading barcode file.")
# check fastq files
fqfiles = sorted(bckeys.keys())
fqok = [tagdigger_fun.isFastq(f) for f in fqfiles]
if not all(fqok):
    print("Cannot read the following as FASTQ files:")
    print([fqfiles[i] for i in range(len(fqfiles)) if not fqok[i]])
    raise Exception("Cannot read all FASTQ files.")

# check binary SNPs for numeric genotype output
if args.outputgen != None:
    if set([t[-1] for t in tags[0]]) != {'0', '1'}:
        raise Exception("Cannot output numeric genotypes for non-binary markers.")

# begin tagdigger algorithm
countsdict = dict()
for f in fqfiles:
    countsdict[f] = tagdigger_fun.find_tags_fastq(f, bckeys[f][0], tags[1],
                                                  cutsite=cutsite)
# combine across libraries
combres = tagdigger_fun.combineReadCounts(countsdict, bckeys)

# output
tagdigger_fun.writeCounts(args.outputcounts, combres[1], combres[0], tags[0])
if args.outputgen != None:
    tagdigger_fun.writeDiploidGeno(args.outputgen, combres[1], combres[0], tags[0])
