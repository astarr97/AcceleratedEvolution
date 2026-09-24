import os
import pandas as pd
import sys

# Read arguments
inFile = sys.argv[1]

assert(inFile.endswith(".txt"))

def complement(bp):
    if bp == "A":
        return("T")
    elif bp == "T":
        return("A")
    elif bp == "C":
        return("G")
    else:
        return("C")


out_hum = []
o = open(inFile)
out = open(inFile.replace(".txt", ".bed"), 'w')
for line in o:
    row = line.replace("\n", "").split("\t")
    new_row_hum = [row[0], str(int(row[1])-1), row[1], row[3]+"|"+row[4]]
    out.write("\t".join(new_row_hum) + "\n")
