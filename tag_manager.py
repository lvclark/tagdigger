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
    print("\nNew tags to add to marker database:")
    tags = tagdigger_fun.readTags_interactive()

    SNPdb = None  # database
    while SNPdb == None:
        SNPdb = tagdigger_fun.readMarkerDatabase(input("Name of CSV file containing marker database: ").strip())

    print("Comparing tags...")
    compareDict = tagdigger_fun.compareTagSets(SNPdb[0], tags)

    # optionally include original marker names
    inclOrig = ""
    origColName = ""
    while inclOrig not in {'Y', 'N'}:
        inclOrig = input("Include column containing original marker names? (y/n) ").strip().upper()
    if inclOrig == 'Y':
        while origColName == "":
            origColName = input("Column header for original marker names: ").strip()

    print("\nCounting markers...")
    # get list of existing markers that had a match
    matchedold = sorted([i for i in compareDict.values() if i != None])
    # get list of all old markers
    allold = sorted(SNPdb[1][1].keys())
    # get a list of all new markers (original names; generate new names below)
    allnew = tagdigger_fun.extractMarkers(tags[0])[0]
    # total number of markers to be output
    nMrkr = len(allold) - len(matchedold) + len(allnew)
    minDig = math.ceil(math.log10(nMrkr)) # minimum number of digits
    # get old marker prefix, number of digits, and last number
    lastold = allold[-1]
    numDig = 0
    for i in [j*-1 for j in range(1, len(lastold))]:
        if lastold[i] in set('0123456789'):
            numDig += 1
        else:
            break
    Prefix = lastold[:-numDig]
    highestNum = int(lastold[-numDig:])
    startingNum = highestNum + 1

    # ask user about keeping the previous numbering system
    print('Last marker name in existing database is {}.'.format(lastold))
    print('Prefix is {}, number of digits is {}, and new markers will be numbered starting {}.'.format(Prefix, numDig, startingNum))
    prefixChoice = input('\nPress enter to keep the prefix {}, or type different prefix to use with new markers: '.format(Prefix)).strip()
    if prefixChoice != '':
        Prefix = prefixChoice

    print('\nTotal number of markers is {}.'.format(nMrkr))
    print('Minimum number of digits is {}.'.format(minDig))
    digChoice = 'a'
    while (not set(digChoice) < set('0123456789')) or (numDig < minDig):
        digChoice = input('\nPress enter to keep {} as the number of digits, or enter a new number: '.format(numDig)).strip()
        if set(digChoice) < set('0123456789') and len(digChoice) > 0:
            numDig = int(digChoice)

    numChoice = 'a'
    while (not set(numChoice) < set('0123456789')) or \
          ("{}{:0{width}}".format(Prefix, startingNum, width = numDig) in allold): # make sure marker names won't be duplicated
        numChoice = input('\nPress enter to start numbering from {}, or enter a different starting number: '.format(startingNum)).strip()
        if set(numChoice) < set('0123456789') and len(numChoice) > 0:
            startingNum = int(numChoice)

    # generate new names for unmatched markers
    print("\nGenerating new marker names...")
    currentNum = startingNum
    unmatchednew = [] # to contain new names of new markers
    for m in allnew:
        if compareDict[m] == None:
            thisnewname = "{}{:0{width}}".format(Prefix, currentNum, width = numDig)
            currentNum += 1
            unmatchednew.append(thisnewname)
            compareDict[m] = thisnewname
    print('{} out of {} markers are new.'.format(len(unmatchednew), len(allnew)))

    # add sequences to tag database
    print("Adding new sequences to tag database...")
    tagsNEW = [[], []] # later combine these two lists with those in SNPdb[0]
    for t in range(len(tags[0])):
        thistagname = tags[0][t]
        thismarker = thistagname[:thistagname.find('_')]
        thismarkerNEW = compareDict[thismarker]
        if thismarkerNEW in unmatchednew:
            thisallelePlusUnd = thistagname[thistagname.rfind('_'):]
            tagsNEW[0].append(thismarkerNEW + thisallelePlusUnd) # new tag name
            tagsNEW[1].append(tags[1][t]) # tag sequence

    # optionally export FASTA of new sequences only
    exptFA = ""
    while exptFA not in {'Y', 'N'}:
        exptFA = input("Make FASTA file of new tags, to use with alignment software? (y/n): ").strip().upper()
    if exptFA == 'Y':
        FAfile = ''
        while FAfile == '':
            FAfile = input("Name for FASTA file: ").strip()
        tagdigger_fun.exportFasta(FAfile, tagsNEW[0], tagsNEW[1])

    # option to also add any tabular data
    addTab = ""
    addTable = None
    while addTab not in {'Y', 'N'}:
        addTab = input("\nAdd additional columns to database, referenced by original marker names? (y/n) ").strip().upper()
    if addTab == 'Y':
        while addTable == None:
            addTable = tagdigger_fun.readTabularData(input("Name of CSV file with additional columns: ").strip(),
                                                     markerDict = compareDict)
        newheaders = addTable[0]
        oldheaders = SNPdb[1][0]
        if len(set(newheaders) & set(oldheaders)) > 0:
            conflChoice = ""
            print('What should be done if conflicting data are found?')
            while conflChoice not in {'O', 'N'}:
                conflChoice = input('o = use old values, n = use new values :').strip().upper()
            if conflChoice == 'O':
                combinedTables = tagdigger_fun.consolidateExtraCols([addTable, SNPdb[1]])
            if conflChoice == 'N':
                combinedTables = tagdigger_fun.consolidateExtraCols([SNPdb[1], addTable])
        else:
            combinedTables = [SNPdb[1], addTable]
    else:
        combinedTables = [SNPdb[1]]

    outfile = ''
    while outfile == '':
        outfile = input("\nName of CSV file for marker database output: ").strip()

    # add original marker names
    if inclOrig == 'Y':
        combinedTables.append([[origColName], {compareDict[k]: [k] for k in compareDict.keys()}])

    # output
    print('\nMaking merged tag sequences...')
    mymerged = tagdigger_fun.mergedTagList([SNPdb[0][0] + tagsNEW[0], SNPdb[0][1] + tagsNEW[1]])
    if mymerged == None:
        print("Please check your input and then re-run the program.")
    else:
        print('Writing file...')
        tagdigger_fun.writeMarkerDatabase(outfile,
                                          mymerged[0], mymerged[1], combinedTables)


# Add alignment data to database
if whichprog == '3':
    SNPdb = None  # database
    while SNPdb == None:
        SNPdb = tagdigger_fun.readMarkerDatabase(input("Name of CSV file containing marker database: ").strip())
    # optionally make fasta file for alignment
    exptFA = ""
    while exptFA not in {'Y', 'N'}:
        exptFA = input("\nMake FASTA file of all tags, to use with alignment software? (y/n): ").strip().upper()
    if exptFA == 'Y':
        FAfile = ''
        while FAfile == '':
            FAfile = input("Name for FASTA file: ").strip()
        tagdigger_fun.exportFasta(FAfile, SNPdb[0][0], SNPdb[0][1])

    # set up dictionary for actual SNP positions, if desired
    varInput = ""
    while varInput not in {'Y', 'N'}:
        varInput = input("\nCalculate actual sites of SNPs, in addition to tag alignment position? (y/n): ").strip().upper()
    if varInput == 'Y':
        print("Variable sites will only be output if there is a single variable site per marker.")
        varDict = tagdigger_fun.varSitesByMarker(SNPdb[0][0], SNPdb[0][1])
    else:
        varDict = None

    # read Bowtie2 output
    bt = None
    while bt == None:
        bt = tagdigger_fun.readSAM(input("\nName of SAM file containing alignment data: ").strip(),
                                   varDict = varDict)

    # make column names
    btColNames = ['','','']
    while btColNames[0] == '':
        btColNames[0] = input('\nName for output column containing chromosome names: ').strip()
    while btColNames[1] == '':
        btColNames[1] = input('Name for output column containing alignment positions: ').strip()
    while btColNames[2] == '':
        btColNames[2] = input('Name for output column containing alignment qualities: ').strip()
    if varInput == 'Y':
        btColNames.append('')
        while btColNames[3] == '':
            btColNames[3] = input('Name for output column containing variable site positions: ').strip()

    # clean up variable sites if needed
    if varInput == 'Y':
        btMarkers = bt.keys() # marker names
        btOut = dict()
        for k in btMarkers:
            thisrow = bt[k]
            outrow = list(thisrow[0:3])
            if len(thisrow[3]) == 1:
                outrow.append(thisrow[3][0])
            else:
                outrow.append("")
            btOut[k] = outrow
    else:
        btOut = bt

    outfile = ''
    while outfile == '':
        outfile = input("\nName of CSV file for marker database output: ").strip()

    # output
    print('\nRemaking merged tag sequences...')
    mymerged = tagdigger_fun.mergedTagList(SNPdb[0])
    print('Writing file...')
    tagdigger_fun.writeMarkerDatabase(outfile,
                                      mymerged[0], mymerged[1],
                                      [SNPdb[1], [btColNames,btOut]])

# Start new database
if whichprog == '4':
    markers = None
    while markers == None:
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

    print('\nOptions for exporting SNP database:')
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

input("\nPress enter to quit.")
