import tagdigger_fun
import math
import csv

# Welcome message
print('''
        TagDigger v. 0.0 Tag Manager
         Copyright Lindsay V. Clark
 Released under GNU General Public License v3
''')

# Need to confirm that marker names contain no spaces, for
# alignment purposes.
# Adjust function for finding adapter sequence so that
# it doesn't look for barcoded adapter when looking at tags.

# set working directory for finding files
tagdigger_fun.set_directory_interactive()

print('''
	Options are:
1. Look up markers by sequence in existing database
2. Add markers to existing database
3. Add alignment data to database
4. Start new database
''')

whichprog  = "0"
while whichprog not in set('1234'):
    whichprog = input("Select option: ").strip()

# Look up markers in existing database
if whichprog == '1':
    print("\nTags to look up in marker database:")
    tags = tagdigger_fun.readTags_interactive() # new markers to look up

    SNPdb = None  # database to query
    while SNPdb == None:
        SNPdb = tagdigger_fun.readMarkerDatabase(input("Name of CSV file containing marker database: ").strip())

    print("Comparing tags...")
    compareDict = tagdigger_fun.compareTagSets(SNPdb[0], tags)

    # choice of which extra columns to include in output; generate numerical index
    inclExtr = ''
    while inclExtr not in {'A', 'S', 'N'}:
        inclExtr = input('''Include additional columns from the database in the table?
a = include all, s = select which to include, n = include none: ''').strip().upper()
    if inclExtr == 'A':
        extracol = list(range(len(SNPdb[1][0])))
    if inclExtr == 'N':
        extracol = []
    if inclExtr == 'S':
        extracol = []
        for i in range(len(SNPdb[1][0])):
            inclThis = ''
            while inclThis not in {'Y', 'N'}:
                inclThis = input("Include {}? (y/n) ".format(SNPdb[1][0][i])).strip().upper()
            if inclThis == 'Y':
                extracol.append(i)

    nHeader = len(SNPdb[1][0]) # number of extra cols, including those not to write

    # file to write
    outfile = ''
    while outfile == '':
        outfile = input("File name for CSV output: ").strip()
    with open(outfile, 'w', newline = '') as outcon:
        cw = csv.writer(outcon)
        # write header row
        cw.writerow(['Query', 'Marker name'] + \
                    [SNPdb[1][0][i] for i in range(nHeader) if i in extracol])
        # write marker rows
        for q in sorted(compareDict.keys()): # loop through query marker names
            dbmarker = compareDict[q]
            if dbmarker == None:
                cw.writerow([q, ''] + ['' for i in range(len(extracol))])
            else:
                cw.writerow([q, dbmarker] + \
                            [SNPdb[1][1][dbmarker][i] for i in range(nHeader) if i in extracol])


# Add markers to existing database
if whichprog == '2':
    tags = tagdigger_fun.readTags_interactive()

    SNPdb = None  # database
    while SNPdb == None:
        SNPdb = tagdigger_fun.readMarkerDatabase(input("Name of CSV file containing marker database: ").strip())
# option to also add any tabular data
# option to make fasta file for unaligned  markers
# what to do with conflicting tabular data
    pass

# Add alignment data to database
if whichprog == '3':
    SNPdb = None  # database
    while SNPdb == None:
        SNPdb = tagdigger_fun.readMarkerDatabase(input("Name of CSV file containing marker database: ").strip())
# optionally make fasta file for alignment
# read Bowtie2 output
    pass

# Start new database
if whichprog == '4':
    tags = tagdigger_fun.readTags_interactive()
    print("Creating merged tag strings for markers...\n")
    markers = tagdigger_fun.mergedTagList(tags)
    nMrkr = len(markers[0]) # number of markers
    minDig = math.ceil(math.log10(nMrkr)) # minimum number of digits

    print('''Markers will be given names in the format Abcde000001.''')
    print('''It is recommended that marker names not include spaces.''')
    mrkrPrefix = ""
    while mrkrPrefix == "":
        mrkrPrefix = input('Prefix for marker names to output ("Abcde" in the above example): ').strip()
    numDig = 0
    while numDig < minDig:
        numDig = int(input('Number of digits for numbering markers (6 in the above example): ').strip())

    # generate marker names
    markerNames = ["{}{:0{width}}".format(mrkrPrefix, i, width = numDig) for i in range(1, nMrkr+1)]
    
    # optionally make fasta file for alignment
    exptFA = ""
    while exptFA not in {'Y', 'N'}:
        exptFA = input("Make FASTA file of tags to use with alignment software? (y/n): ").strip().upper()
    if exptFA == 'Y':
        FAfile = ''
        while FAfile == '':
            FAfile = input("Name for FASTA file: ").strip()
        tagdigger_fun.exportFasta2(FAfile, markerNames, markers[1])

    # optionally include original marker names
    inclOrig = ""
    origColName = ""
    while inclOrig not in {'Y', 'N'}:
        inclOrig = input("Include column containing original marker names? (y/n) ").strip().upper()
    if inclOrig == 'Y':
        while origColName == "":
            origColName = input("Column header for original marker names: ").strip()

    # optionally add tabular data
    addTab = ""
    addTable = None
    while addTab not in {'Y', 'N'}:
        addTab = input("Add additional columns of data, referenced by original marker names? (y/n) ").strip().upper()
    if addTab == 'Y':
        while addTable == None:
            addTable = tagdigger_fun.readTabularData(input("Name of CSV file with additional columns: ").strip(),
                                                     markerDict = dict(zip(markers[0], markerNames)))

    # output to CSV file
    extracollist = list()
    if addTab == 'Y':
        extracollist.append(addTable)
    if inclOrig == 'Y':
        extracollist.append([[origColName], dict(zip(markerNames, [[m] for m in markers[0]]))])
    outfile = ''
    while outfile == '':
        outfile = input("Name of CSV file for marker database output: ").strip()
    tagdigger_fun.writeMarkerDatabase(outfile,
                                      markerNames, markers[1], extracollist)

input("Press enter to quit.")
