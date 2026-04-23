#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

spec_mat="VIP/SpeciesMatrix_VIP_abs_logfc_NC_FiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90"
back_mat="VIP/BackgroundMatrix_FiltWGS_VIP_NC_PhyloP1_SpecSup250_NCM_Per90"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 1 250 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_FiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90_PerGene VIP abs_logfc $spec_mat $back_mat 0,1000

spec_mat="VIP/SpeciesMatrix_VIP_abs_logfc_NC_NoFiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90"
back_mat="VIP/BackgroundMatrix_NoFiltWGS_VIP_NC_PhyloP1_SpecSup250_NCM_Per90"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 0 1 250 0 1 0 90 NC Outputs/VIP_AbsLogfc_NC_NoFiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90_PerGene VIP abs_logfc $spec_mat $back_mat 0,1000
