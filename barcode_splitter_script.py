#!/usr/bin/python3

# non-interactive script for doing barcode splitting

import tagdigger_fun
import argparse

parser = argparse.ArgumentParser(description="TagDigger v. 1.1 barcode splitter command line script by Lindsay V. Clark")
parser.add_argument('-b', '--barcodefile', help = 'Name of barcode key file',
                    required = True)
parser.add_argument('-a', '--adapter', help = 'Name of the adapter set', required = True,
                    choices = sorted(tagdigger_fun.adapters.keys()))

args = parser.parse_args()

bckeys = tagdigger_fun.readBarcodeKeyfile(args.barcodefile, forSplitter = True)

adapter = tagdigger_fun.adapters[args.adapter]

enzyme = args.adapter[:args.adapter.find("-")]

cutsite = tagdigger_fun.enzymes[enzyme]

fqfiles = sorted(bckeys.keys())

for f in fqfiles:
  tagdigger_fun.barcodeSplitter(f, bckeys[f][0], bckeys[f][1], 
                                cutsite = cutsite, adapter = adapter)
