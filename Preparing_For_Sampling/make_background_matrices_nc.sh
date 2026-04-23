#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=128GB

python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 1 -100 0 0 Shared/BackgroundMatrix_FiltWGS_PhyloP_NC_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 0 -100 0 0 Shared/BackgroundMatrix_NoFiltWGS_PhyloP_NC_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 1 -100 250 0 Shared/BackgroundMatrix_FiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 0 -100 250 0 Shared/BackgroundMatrix_NoFiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_NCM_Per0 HumChp_AccelEvolInput.Final.txt

python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 1 -100 0 1 PhyloP/BackgroundMatrix_FiltWGS_PhyloP_NC_PhyloP-100_SpecSup0_YCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 0 -100 0 1 PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_NC_PhyloP-100_SpecSup0_YCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 1 -100 250 1 PhyloP/BackgroundMatrix_FiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_YCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.NC.Final.sort.txt NC 0 -100 250 1 PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_NC_PhyloP-100_SpecSup250_YCM_Per0 HumChp_AccelEvolInput.Final.txt
