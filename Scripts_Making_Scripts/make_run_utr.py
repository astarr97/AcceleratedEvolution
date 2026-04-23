sim_num = [100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000]

def write_header(out):
    out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hbfraser\n#SBATCH --mem=16GB\n\n")

for i in sim_num:
    #Write out for the 3' UTR stuff
    #out = open("run_3utr_paramset1_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Filter WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    #out.close()
    
    #out = open("run_3utr_paramset2_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Not filtering WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    #out.close()
    
    out = open("run_3utr_paramset3_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    out.close()

    out = open("run_3utr_paramset4_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Not filtering WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    out.close()
    
    #Write out for the 5' UTR stuff
    #out = open("run_5utr_paramset1_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Filter WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    #out.close()
    
    #out = open("run_5utr_paramset2_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Not filtering WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 0 5UTR Outputs/PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    #out.close()
    
    out = open("run_5utr_paramset3_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    out.close()

    out = open("run_5utr_paramset4_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Not filtering WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    out.close()

sim_num = list(range(10000, 1010000, 10000))
for i in sim_num:
    #Write out the 3' UTR stuff for HPO
    #out = open("run_3utr_paramset1_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Filter WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    #out = open("run_3utr_paramset2_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Not filtering WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    out = open("run_3utr_paramset3_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()

    out = open("run_3utr_paramset4_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Not filtering WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()
    
    
    
    #Write out the 3' UTR stuff for GOBP
    #out = open("run_3utr_paramset1_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Filter WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    #out = open("run_3utr_paramset2_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Not filtering WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    out = open("run_3utr_paramset3_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()

    out = open("run_3utr_paramset4_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Not filtering WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()
    
    
    #Write out the 5' UTR stuff for HPO
    #out = open("run_5utr_paramset1_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Filter WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    #out = open("run_5utr_paramset2_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Not filtering WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 0 5UTR Outputs/PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    out = open("run_5utr_paramset3_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()

    out = open("run_5utr_paramset4_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Not filtering WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()

    
    #Write out the 5' UTR stuff for GOBP
    #out = open("run_5utr_paramset1_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Filter WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    #out = open("run_5utr_paramset2_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    #write_header(out)
    #out.write('#Not filtering WGS, no other parameters\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0"\n')
    #out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 0 5UTR Outputs/PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    #out.close()
    
    out = open("run_5utr_paramset3_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()

    out = open("run_5utr_paramset4_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Not filtering WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()
    
"""for i in sim_num:
    out = open("run_3utr_paramset3_" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    out.close()

    out = open("c" + str(i - 100000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 100000), str(i)]) + "\n")
    out.close()

sim_num = list(range(10000, 1010000, 10000))
for i in sim_num:
    out = open("run_3utr_paramset3_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()
    out = open("run_3utr_paramset3_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()
    out = open("run_5utr_paramset3_hpo_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()
    out = open("run_5utr_paramset3_gobp_" + str(i - 10000) + "-" + str(i) + ".sh", 'w')
    write_header(out)
    out.write('#Filter WGS, SpecSup250\nspec_mat="PhyloP/SpeciesMatrix_PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"\nback_mat="PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0"\n')
    out.write("python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 0 0 0 5UTR Outputs/PhyloP_5UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat " + ",".join([str(i - 10000), str(i)]) + " /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100\n")
    out.close()"""
    