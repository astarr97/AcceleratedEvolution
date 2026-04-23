import os
import pandas as pd
import sys

# Read arguments
inFile = sys.argv[1]
outFile = sys.argv[2]

def complement(bp):
    if bp == "A":
        return("T")
    elif bp == "T":
        return("A")
    elif bp == "C":
        return("G")
    else:
        return("C")


o = open(inFile)
out = open(outFile, 'w')
for line in o:
    row = line.replace("\n", "").split("\t")
    new_row_hum = [row[0], str(int(row[1])-1), str(row[1]), row[3]+"|"+row[4]]
    out.write("\t".join(new_row_hum) + "\n")
o.close()
out.close()
