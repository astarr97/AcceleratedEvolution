o = open("HumChp_AccelEvolInput.ToFixGene.bed")
out = open("HumChp_AccelEvolInput.Final.txt", 'w')
out.write("Position\tHum|Chp|Gor\tDerived\tPhyloP447\tSpecSup447\tPhastCons447\tCategory\tKeptAfterFilt\tAncTrinuc\tDerTrinuc\tMutCat\tNearestGene\tNearestDist\n")

for line in o:
    l = line.split("\t")
    if l[8] == "NC" or l[8] == "NC_Syn" or l[8] == ".":
        if "-" in l[13]:
            if l[13] == "10-Mar":
                l[13] = "MARCHF10"
            elif l[13] == "11-Mar":
                l[13] = "MARCHF11"
            elif l[13] == "1-Mar":
                if l[0] == "chr4":
                    l[13] = "MARCHF1"
                elif l[0] == "chr1":
                    l[13] = "MTARC1"
            elif l[13] == "2-Mar":
                l[13] = "MARCHF2"
            elif l[13] == "3-Mar":
                l[13] = "MARCHF3"
            elif l[13] == "4-Mar":
                l[13] = "MARCHF4"
            elif l[13] == "5-Mar":
                l[13] = "MARCHF5"
            elif l[13] == "6-Mar":
                l[13] = "MARCHF6"
            elif l[13] == "7-Mar":
                l[13] = "MARCHF7"
            elif l[13] == "8-Mar":
                l[13] = "MARCHF8"
            elif l[13] == "9-Mar":
                l[13] = "MARCHF9"
            if l[13] == "10-Sep":
                l[13] = "SEPTIN10"
            elif l[13] == "11-Sep":
                l[13] = "SEPTIN11"
            elif l[13] == "1-Sep":
                l[13] = "SEPTIN1"
            elif l[13] == "2-Sep":
                l[13] = "SEPTIN2"
            elif l[13] == "3-Sep":
                l[13] = "SEPTIN3"
            elif l[13] == "4-Sep":
                l[13] = "SEPTIN4"
            elif l[13] == "5-Sep":
                l[13] = "SEPTIN5"
            elif l[13] == "6-Sep":
                l[13] = "SEPTIN6"
            elif l[13] == "7-Sep":
                l[13] = "SEPTIN7"
            elif l[13] == "8-Sep":
                l[13] = "SEPTIN8"
            elif l[13] == "9-Sep":
                l[13] = "SEPTIN9"
            elif l[13] == "12-Sep":
                l[13] = "SEPTIN12"
            elif l[13] == "13-Sep":
                l[13] = "SEPTIN13"
            elif l[13] == "14-Sep":
                l[13] = "SEPTIN14"
            elif l[13] == "15-Sep":
                l[13] = "SEPTIN15"
            elif l[13] == "1-Dec":
                l[13] = "DELEC1"
    else:
        if (l[8] == "5UTR_Syn" or l[8] == "5UTR") and l[-2] != "5UTR":
            l[8] = "."
        if (l[8] == "3UTR_Syn" or l[8] == "3UTR") and l[-2] != "3UTR":
            l[8] = "."
        l[13] = l[18]
    out.write("\t".join([l[0] + ":" + l[2]] + l[3:15]) + "\n")
