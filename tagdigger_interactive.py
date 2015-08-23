# script for running an interactive session with TagDigger.

import tagdigger_fun
import os

# Welcome message
print('''
                  TagDigger v. 0.0
             Copyright Lindsay V. Clark
    Released under GNU General Public License v3
    ''')

# Choose enzyme
knownenzymes = sorted(tagdigger_fun.enzymes.keys())
print("Known restriction enzymes are:")
enzlines = ""
for i in range(len(knownenzymes)):
    enzlines += knownenzymes[i]
    if i % 8 == 7:
        enzlines += "\n"
    else:
        enzlines += " "
print(enzlines)

print('''
What restriction cut site should be found immediately
after the barcode sequence?  Type the name of one of the
above enzymes, OR type the restriction cut site as it
should appear in the sequence data (i.e. not including
bases before the beginning of the overhang) using
characters ACGTRYSWKMBDHVN (IUPAC codes for ambiguous
nucleotides).
''')
enzdone = False
while not enzdone:
    enzchoice = input("Restriction site: ")
    if enzchoice in knownenzymes:
        cutsite = tagdigger_fun.enzymes[enzchoice]
        enzdone = True
    elif set(enzchoice.upper()) <= set('ACGTRYSWKMBDHVN'):
        cutsite = enzchoice.upper()
        enzdone = True
print("Cut site: " + cutsite)

# set working directory for finding files
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
        toKeep = tagdigger_fun.readMarkerNames(input("File name: ").strip())
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
''')
tagfunctions = {'1': tagdigger_fun.readTags_UNEAK_FASTA,
                '2': tagdigger_fun.readTags_Merged,
                '3': tagdigger_fun.readTags_Columns,
                '4': tagdigger_fun.readTags_Rows}

# choose format and read tag file
tags = None
while tags == None:
    thischoice = '0'
    while thischoice not in {'1', '2', '3', '4'}:
        thischoice = input("Enter the number of the format of your tag file: ").strip()
    tagfile = input("Enter the file name: ").strip()
    tags = tagfunctions[thischoice](tagfile, toKeep = toKeep)
    print('')
# summarize results
print("{} tag sequences read.\n".format(len(tags[1])))

# sanitize
tags = tagdigger_fun.sanitizeTags(tags)
print("{} tag sequences remain.\n".format(len(tags[1])))

# read in key file
bckeys = None
while bckeys == None:
    bckeys = tagdigger_fun.readBarcodeKeyfile(input("Name of key file with barcodes: ").strip())
# summarize
fqfiles = sorted(bckeys.keys())
for f in fqfiles:
    print("File {}: {} barcodes".format(f, len(bckeys[f][0])))
print("")

# check that FASTQ files are okay
fqok = [tagdigger_fun.isFastq(f) for f in fqfiles]
while not all(fqok):
    print("Cannot read the following as FASTQ files:")
    for f in range(len(fqok)):
        if not fqok[f]:
            print(fqfiles[f])
    thischoice = '0'
    while thischoice not in {'1', '2', '3'}:
        thischoice = input('''
Press 1 to re-read key file, 2 to search for FASTQ files in a different
directory, or 3 to try reading the same FASTQ files again: ''').strip()
    if thischoice == '1':
        bckeys = None
        while bckeys == None:
            bckeys = tagdigger_fun.readBarcodeKeyfile(input("\nName of key file with barcodes: "))
        fqfiles = sorted(bckeys.keys())
        for f in fqfiles:
            print("File {}: {} barcodes".format(f, len(bckeys[f])))
        print("")
    if thischoice == '2':
        dirchoice = ""
        while not os.path.isdir(dirchoice):
            dirchoice = input("New directory: ")
        os.chdir(dirchoice)
    fqok = [tagdigger_fun.isFastq(f) for f in fqfiles]

# Get file name for output
countsfile = ""
while countsfile =="":
    countsfile = input("\nFile name for output of read counts: ").strip()

# check whether these are binary SNPs, and ask whether to output as numeric genotypes
genofile = ""
if set([t[-1] for t in tags[0]]) == {'0', '1'}:
    thischoice = ""
    while thischoice not in {'Y', 'N'}:
        thischoice = input("\nOutput CSV of diploid numeric genotypes? Y/N ").strip().upper()
    if thischoice == 'Y':
        while genofile == "":
            genofile = input("File name for output of genotypes: ").strip()

# Run the tag search
input("\nPress enter to begin processing FASTQ files.")
countsdict = dict()
for f in fqfiles:
    countsdict[f] = tagdigger_fun.find_tags_fastq(f, bckeys[f][0], tags[1],
                                                  cutsite=cutsite)

# combine across libraries
combres = tagdigger_fun.combineReadCounts(countsdict, bckeys)

# Confirm directory for output

# Output tag counts if desired
tagdigger_fun.writeCounts(countsfile, combres[1], combres[0], tags[0])

# Output diploid genotypes if desired
if genofile != "":
    tagdigger_fun.writeDiploidGeno(genofile, combres[1], combres[0], tags[0])

input("\nPress enter to quit.")
    
