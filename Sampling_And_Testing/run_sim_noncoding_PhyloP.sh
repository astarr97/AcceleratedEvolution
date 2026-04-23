#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=128GB

#Now with WGS filtering
spec_mat="PhyloP/SpeciesMatrix_PhyloP_NC_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_FiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_NCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100

spec_mat="PhyloP/SpeciesMatrix_PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_NoFiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_NCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 0 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100

spec_mat="PhyloP/SpeciesMatrix_PhyloP_NC_FiltWGS_PhyloP-100_SpecSup250_YCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_FiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_YCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_YCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_YCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_YCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100

spec_mat="PhyloP/SpeciesMatrix_PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_YCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_NoFiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_YCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 0 0 NC Outputs/PhyloP_NC_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100
