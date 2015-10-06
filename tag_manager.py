import tagdigger_fun
import math

# Welcome message
print('''
TagDigger v. 0.0 Tag Manager
Copyright Lindsay V. Clark
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
    tags = tagdigger_fun.readTags_interactive()
# Export table with marker names from input in first column,
# user-selected columns from database in remaining columns.

# Add markers to existing database
if whichprog == '2':
    tags = tagdigger_fun.readTags_interactive()
# option to also add any tabular data
# option to make fasta file for unaligned  markers
# what to do with conflicting tabular data
    pass

# Add alignment data to database
if whichprog == '3':
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
