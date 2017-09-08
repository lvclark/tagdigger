import tagdigger_fun
import os

# Welcome message
print('''
     TagDigger v. 1.1 Barcode Splitter
        Copyright Lindsay V. Clark
Released under GNU General Public License v3
''')

# Choose enzyme
knownenzymes = ["NsiI", "PstI"]
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
characters ACGT.
''')
enzdone = False
while not enzdone:
    enzchoice = input("Restriction site: ")
    if enzchoice in knownenzymes:
        cutsite = tagdigger_fun.enzymes[enzchoice]
        enzdone = True
    elif set(enzchoice.upper()) <= set('ACGT'):
        cutsite = enzchoice.upper()
        enzdone = True
print("Cut site: " + cutsite)

# choose adapter set
print("\nKnown adapter sets:")
adaptersets = sorted(tagdigger_fun.adapters.keys())
for a in adaptersets:
    if enzchoice not in knownenzymes or enzchoice in a:
        print(a)
print("")
adaptchoice = ""
while adaptchoice not in adaptersets:
    adaptchoice = input("Choose an adapter set: ").strip()

# set working directory for finding files
tagdigger_fun.set_directory_interactive()

# read in key file
bckeys = None
while bckeys == None:
    bckeys = tagdigger_fun.readBarcodeKeyfile(input("\nName of key file with barcodes: ").strip(), forSplitter=True)
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

# ask if MD5 checksums should be calculated
print('')
md5choice = ""
while md5choice not in {'Y', 'N'}:
    md5choice = input("Create a CSV file of MD5 checksums? (y/n) ").strip().upper()
if md5choice == 'Y':
    md5outfile = ''
    while md5outfile == '':
        md5outfile = input("Name of output CSV file to contain MD5 checksums: ").strip()
        

input("\nPress enter to begin processing files.")

# go through files
for f in fqfiles:
    tagdigger_fun.barcodeSplitter(f, bckeys[f][0], bckeys[f][1], cutsite = cutsite,
                                  adapter = tagdigger_fun.adapters[adaptchoice])

# calculate MD5 checksums
if md5choice == 'Y':
    filelist = []
    for f in fqfiles:
        filelist += bckeys[f][1]
    tagdigger_fun.writeMD5sums(filelist, md5outfile)
    

input("\nPress enter to quit.")
