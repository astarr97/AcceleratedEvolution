#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=128GB

python make_background_matrices.py Background_AccelEvolInput.5UTR.Final.sort.txt 5UTR 1 -100 0 0 PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.5UTR.Final.sort.txt 5UTR 0 -100 0 0 PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.5UTR.Final.sort.txt 5UTR 1 -100 250 0 PhyloP/BackgroundMatrix_FiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.5UTR.Final.sort.txt 5UTR 0 -100 250 0 PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_5UTR_PhyloP-100_SpecSup250_NCM_Per0 HumChp_AccelEvolInput.Final.txt

python make_background_matrices.py Background_AccelEvolInput.3UTR.Final.sort.txt 3UTR 1 -100 0 0 PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.3UTR.Final.sort.txt 3UTR 0 -100 0 0 PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.3UTR.Final.sort.txt 3UTR 1 -100 250 0 PhyloP/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0 HumChp_AccelEvolInput.Final.txt
python make_background_matrices.py Background_AccelEvolInput.3UTR.Final.sort.txt 3UTR 0 -100 250 0 PhyloP/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup250_NCM_Per0 HumChp_AccelEvolInput.Final.txt
