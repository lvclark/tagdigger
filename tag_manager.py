import tagdigger_fun

# Welcome message
print('''
TagDigger v. 0.0 Tag Manager
Copyright Lindsay V. Clark
''')

# Need function to export tags in merged format.
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
# Export table with marker names from input in first column,
# user-selected columns from database in remaining columns.
    pass

# Add markers to existing database
if whichprog == '2':
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
# optionally make fasta file for alignment 
# optionally add tabular data.
    pass

input("Press enter to quit.")
