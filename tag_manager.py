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
    markerNames = ["{}{:0>}".format(mrkrPrefix, i) for i in range(1, nMrkr+1)]
    
# optionally make fasta file for alignment 
# optionally add tabular data.

input("Press enter to quit.")
