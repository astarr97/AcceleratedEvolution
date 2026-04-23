#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH -p hbfraser
#SBATCH --mem=64GB

#Name of the .txt file (always the same)
#file = sys.argv[1]

#Whether to filter out the sites that might be polymorphic (1 is yes, 0 is no)
#filt_wgs = int(sys.argv[2])

#Variant effect prediction effect size cutoff
#p_cut = float(sys.argv[3])

#Species support cutoff
#spec_sup = float(sys.argv[4])

#Whether to control for the type of mutation (1 is yes, 0 is no)
#cont_mut = int(sys.argv[5])

#Whether to bin by variant effect prediction when doing the assignment
#cont_var = int(sys.argv[6])

#Whether we are controlling for the predictions in the background
#back_cont_var = int(sys.argv[7])

#Percentile cutoff to consider an element "accessible"
#percent_cut = float(sys.argv[8])

#Must be NC, 3UTR, 3UTR, or Mis
#variant_cat = sys.argv[9]

#Folder to write the matrices out to
#out_folder = sys.argv[10]

#Either the name of a deep learning prediction set or "PhyloP447"
#var_pred = sys.argv[11]

#Variant effect prediction metric to use, can be "abs logfc" or "logfc" if it is deep learning, can only be "PhyloP447" otherwise
#metric = sys.argv[12].replace("_", " ")

#Folder where the species probabilities are stored
#spec_prob_folder = sys.argv[13]

#Folder where the background probabilities are stored
#back_prob_folder = sys.argv[14]

#Whether we are doing this at the gene set level
#Input should be of the form file,maximum number of genes per category,minimum number of genes per category
#try:
#    gene_set_string = sys.argv[15]
#    gene_set = gene_set_string.split(",")[0]
#    min_genes = int(gene_set_string.split(",")[1])
#    max_genes = int(gene_set_string.split(",")[2])
#except:
#    gene_set = 0

#Filter WGS, no other parameters
spec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100

spec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100

spec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_FiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_FiltWGS_PhyloP-100_SpecSup0_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100

spec_mat="PhyloP/SpeciesMatrix_PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0"
back_mat="Shared/BackgroundMatrix_NoFiltWGS_PhyloP_3UTR_PhyloP-100_SpecSup0_NCM_Per0"
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene PhyloP447 PhyloP447 $spec_mat $back_mat 0,1000000
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_HPO PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/HPO_AccelEvol_Input_FiltForAccelEvol.txt,15,100
python simulate_new_humchp.py HumChp_AccelEvolInput.Final.txt 1 -100 0 0 0 0 0 3UTR Outputs/PhyloP_3UTR_NoFiltWGS_PhyloP-100_SpecSup250_NCM_NCV_Per0_PerGene_GOBP PhyloP447 PhyloP447 $spec_mat $back_mat 0,10000 /home/groups/hbfraser/astarr_scripts/AccelConv/GeneSets/GOBP_AccelEvol_Input_FiltForAccelEvol.txt,15,100
