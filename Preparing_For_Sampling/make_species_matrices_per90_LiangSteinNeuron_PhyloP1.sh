#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

python make_species_matrices.py HumChp_AccelEvolInput.Final.txt 1 -100 250 0 1 90 NC LiangSteinNeuron/SpeciesMatrix_LiangSteinNeuron_abs_logfc_NC_FiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90 LiangSteinNeuron abs_logfc
python make_species_matrices.py HumChp_AccelEvolInput.Final.txt 0 -100 250 0 1 90 NC LiangSteinNeuron/SpeciesMatrix_LiangSteinNeuron_abs_logfc_NC_NoFiltWGS_PhyloP1_SpecSup250_YCM_NCV_Per90 LiangSteinNeuron abs_logfc
