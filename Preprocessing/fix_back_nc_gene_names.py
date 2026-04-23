o = open("Background_AccelEvolInput.NC.Final.sort.temp.txt")
out = open("Background_AccelEvolInput.NC.Final.sort.FixedGeneNames.txt", 'w')

for line in o:
    l = line.split("\t")
    if "-" in l[3]:
        if l[3] == "10-Mar":
            l[3] = "MARCHF10"
        elif l[3] == "11-Mar":
            l[3] = "MARCHF11"
        elif l[3] == "1-Mar":
            if l[0].split(":")[0] == "chr4":
                l[3] = "MARCHF1"
            elif l[0].split(":")[0] == "chr1":
                l[3] = "MTARC1"
        elif l[3] == "2-Mar":
            l[3] = "MARCHF2"
        elif l[3] == "3-Mar":
            l[3] = "MARCHF3"
        elif l[3] == "4-Mar":
            l[3] = "MARCHF4"
        elif l[3] == "5-Mar":
            l[3] = "MARCHF5"
        elif l[3] == "6-Mar":
            l[3] = "MARCHF6"
        elif l[3] == "7-Mar":
            l[3] = "MARCHF7"
        elif l[3] == "8-Mar":
            l[3] = "MARCHF8"
        elif l[3] == "9-Mar":
            l[3] = "MARCHF9"
        if l[3] == "10-Sep":
            l[3] = "SEPTIN10"
        elif l[3] == "11-Sep":
            l[3] = "SEPTIN11"
        elif l[3] == "1-Sep":
            l[3] = "SEPTIN1"
        elif l[3] == "2-Sep":
            l[3] = "SEPTIN2"
        elif l[3] == "3-Sep":
            l[3] = "SEPTIN3"
        elif l[3] == "4-Sep":
            l[3] = "SEPTIN4"
        elif l[3] == "5-Sep":
            l[3] = "SEPTIN5"
        elif l[3] == "6-Sep":
            l[3] = "SEPTIN6"
        elif l[3] == "7-Sep":
            l[3] = "SEPTIN7"
        elif l[3] == "8-Sep":
            l[3] = "SEPTIN8"
        elif l[3] == "9-Sep":
            l[3] = "SEPTIN9"
        elif l[3] == "12-Sep":
            l[3] = "SEPTIN12"
        elif l[3] == "13-Sep":
            l[3] = "SEPTIN13"
        elif l[3] == "14-Sep":
            l[3] = "SEPTIN14"
        elif l[3] == "15-Sep":
            l[3] = "SEPTIN15"
        elif l[3] == "1-Dec":
            l[3] = "DELEC1"
    out.write("\t".join(l))
o.close()
out.close()
