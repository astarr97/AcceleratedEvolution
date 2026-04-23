import os

"""template = open("run_sim_LiangSteinNeuron_test.sh")
c = 1
out = open("run_sim_LiangSteinNeuron_Params" + str(c) + ".sh", 'w')
out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
for line in template:
    if line == "\n":
        out.close()
        c+=1
        out = open("run_sim_LiangSteinNeuron_Params" + str(c) + ".sh", 'w')
        if c == 5 or c == 17 or c == 21:
            out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
        else:
            out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
    else:
        out.write(line)
out.close()

template = open("run_sim_LiangSteinNeuron_test2.sh")
c = 1
out = open("run_sim_LiangSteinNeuron_Params" + str(c) + ".sh", 'w')
out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
for line in template:
    if line == "\n":
        out.close()
        c+=1
        out = open("run_sim_LiangSteinNeuron_YCM_Params" + str(c) + ".sh", 'w')
        if c == 5 or c == 17 or c == 21:
            out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
        else:
            out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
    else:
        out.write(line)
out.close()

c = 1
o = open("Config_ToRun.txt")
driver = open("driver_per90_AbsLogfc.sh", 'w')
for line in o:
    ct = line.replace("\n", "")
    
    template = open("run_sim_Per90_AbsLogfc_template.sh")

    out = open("run_sim_" + ct + "_AbsLogfc_Params_AllPer90.sh", 'w')
    driver.write("sbatch -p hbfraser " + "run_sim_" + ct + "_AbsLogfc_Params_AllPer90.sh\n")
    out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
    for line in template:
        out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
o.close()
driver.close()

c = 1
o = open("Config_ToRun.txt")
driver = open("driver_per90_Logfc.sh", 'w')
for line in o:
    ct = line.replace("\n", "")
    
    template = open("run_sim_Per90_Logfc_template.sh")

    out = open("run_sim_" + ct + "_Logfc_Params_AllPer90.sh", 'w')
    driver.write("sbatch -p hbfraser " + "run_sim_" + ct + "_AbsLogfc_Params_AllPer90.sh\n")
    out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
    for line in template:
        out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
o.close()
driver.close()"""

c = 1
o = open("Config_ToRun.txt")
driver = open("driver_per90_AbsLogfc_PhyloP1.sh", 'w')
for line in o:
    ct = line.replace("\n", "")
    
    template = open("run_sim_Per90_AbsLogfc_template_PhyloP1.sh")

    out = open("run_sim_" + ct + "_AbsLogfc_Params_AllPer90_PhyloP1.sh", 'w')
    driver.write("sbatch -p hbfraser " + "run_sim_" + ct + "_AbsLogfc_Params_AllPer90_PhyloP1.sh\n")
    out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=64GB\n\n")
    for line in template:
        out.write(line.replace("LiangSteinNeuron", ct))
    out.close()
o.close()
driver.close()