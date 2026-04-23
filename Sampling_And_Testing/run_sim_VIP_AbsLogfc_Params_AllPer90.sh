#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

#spec_mat="VIP/SpeciesMatrix_VIP_abs_logfc_NC_FiltWGS_abs_logfc0_SpecSup0_NCM_NCV_Per90"
#back_mat="VIP/BackgroundMatrix_FiltWGS_VIP_NC_PhyloP-100_SpecSup0_NCM_Per90"
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90_PerGene VIP abs_logfc $spec_mat $back_mat 0,1000
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90_PerGene_HPO VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90_PerGene_GOBP VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100

#spec_mat="VIP/SpeciesMatrix_VIP_abs_logfc_NC_NoFiltWGS_abs_logfc0_SpecSup0_NCM_NCV_Per90"
#back_mat="VIP/BackgroundMatrix_NoFiltWGS_VIP_NC_PhyloP-100_SpecSup0_NCM_Per90"
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90_PerGene VIP abs_logfc $spec_mat $back_mat 0,1000
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90_PerGene_HPO VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 0 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per90_PerGene_GOBP VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100


#spec_mat="VIP/SpeciesMatrix_VIP_abs_logfc_NC_FiltWGS_abs_logfc0_SpecSup0_YCM_NCV_Per90"
#back_mat="VIP/BackgroundMatrix_FiltWGS_VIP_NC_PhyloP-100_SpecSup0_NCM_Per90"
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene VIP abs_logfc $spec_mat $back_mat 0,1000
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_HPO VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_GOBP VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100

spec_mat="VIP/SpeciesMatrix_VIP_abs_logfc_NC_NoFiltWGS_abs_logfc0_SpecSup0_YCM_NCV_Per90"
back_mat="VIP/BackgroundMatrix_NoFiltWGS_VIP_NC_PhyloP-100_SpecSup0_NCM_Per90"
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene VIP abs_logfc $spec_mat $back_mat 0,1000
#python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_HPO VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 -100 0 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP-100_SpecSup0_YCM_NCV_Per90_PerGene_GOBP VIP abs_logfc $spec_mat $back_mat 0,1000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol_ForML.txt,15,100
