import os

"""driver = open("matrix_making_driver.sh", 'w')

o = open("Config_ToRun.txt")
for ct in o:
    ct = ct.replace("\n", "")
    if ct not in os.listdir():
        os.mkdir(ct)
    out = open("matrix_making_" + ct + ".sh", 'w')
    template = open("matrix_making_dl.sh")
    for line in template:
        if line.startswith("#") or line == "\n":
            out.write(line)
        else:
            out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
    template.close()
    driver.write("sbatch -p hbfraser " + "matrix_making_" + ct + ".sh" + "\n")
o.close()

o = open("Config_ToRun.txt")
for ct in o:
    ct = ct.replace("\n", "")
    if ct not in os.listdir():
        os.mkdir(ct)
    out = open("matrix_making2_" + ct + ".sh", 'w')
    template = open("matrix_making_dl2.sh")
    for line in template:
        if line.startswith("#") or line == "\n":
            out.write(line)
        else:
            out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
    template.close()
    driver.write("sbatch -p hbfraser " + "matrix_making2_" + ct + ".sh" + "\n")
o.close()

driver = open("matrix_making_driver_Per90.sh", 'w')

o = open("Config_ToRun.txt")
for ct in o:
    ct = ct.replace("\n", "")
    if ct not in os.listdir():
        os.mkdir(ct)
    out = open("matrix_making_Per90_" + ct + ".sh", 'w')
    template = open("make_background_matrices_per90_test.sh")
    for line in template:
        if line.startswith("#") or line == "\n":
            out.write(line)
        else:
            out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
    template.close()
    driver.write("sbatch -p hbfraser " + "matrix_making_Per90_" + ct + ".sh" + "\n")
o.close()"""

driver = open("matrix_making_driver_Per90_PhyloP1.sh", 'w')

o = open("Config_ToRun.txt")
for ct in o:
    ct = ct.replace("\n", "")
    if ct not in os.listdir():
        os.mkdir(ct)
    out = open("matrix_making_Per90_PhyloP1_" + ct + ".sh", 'w')
    template = open("make_background_matrices_per90_PhyloP1_test.sh")
    for line in template:
        if line.startswith("#") or line == "\n":
            out.write(line)
        else:
            out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
    template.close()
    driver.write("sbatch -p hbfraser " + "matrix_making_Per90_PhyloP1_" + ct + ".sh" + "\n")
o.close()
