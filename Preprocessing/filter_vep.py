import sys
import os

for file in os.listdir():
    if ".vcf" in file and "Filt" not in file and "GW_Ens" in file:
        print(file)
        o = open(file)
        out = open(file.replace(".vcf", "_Filt.vcf"), 'x')
        
        for line in o:
            l = line.split("\t")
            if "#" in line:
                out.write(line)
            elif l[2].split("/")[1] == l[4]:
                out.write(line)