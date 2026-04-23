To prepare for the sampling procedure, the first step is to make probability matrices to assign substitutions to a lineage and a gene during sampling
The matrices for lineages are referred to as "species matrices" and the matrices for genes are referred to as "background matrices"
To make the species matrices, we use make_species_matrices.py.  Please see that file for details about implementation and arguments
This needs to be done any time you change a parameter or use a different deep learning model if filtering by accessibility.
In general, we made species matrices for non-coding, 3' UTR, and 5' UTR PhyloP-based analyses with species support > 250, ChromBPNet predictions for 34 cell types, restricting to substitutions in sites in the 90th percentile of accessibility or higher, and then for ChromBPNet predictions again but restricting only to sites with species support > 250 and PhyloP > 1.
If you are controlling for shifts in the PhyloP or predicted absolute log fold-change in CA as we did in the paper, then the output of this is a 32x100 matrix of probabilities.
Each row is one of 32 trinucleotides (32 because AAA is equivalent to TTT etc.) and each column is a PhyloP or predicted absolute log fold-change bin.  The entries are the probability that a substitution would be assigned to the human lineage (it would be 1- this probability for assignment to the chimp lineage).

For background matrices for non-coding or UTR PhyloP, we just use make_background_matrices.py (see script for details).  Examples are shown in make_background_matrices_nc.sh.
The output of this is a 32xN matrix where N is the number of genes and each entry reflects the probability of a substitution in that ancestral trinucleotide being assigned to that gene.

For background matrices for ChromBPNet, we first use prepare_per90_background.py and the code in make_background_matrices_per90_LiangSteinNeuron.sh to get out a set of regions flanking sites with substitutions in the top 90th percentile of accessibility and restrict our background file to only those regions.  We then use make_background_matrices_per90.py, which is very similar but just uses the filtered background input file.
This was done for all 34 cell types separately.

The heads of the input files are given as a guide.
