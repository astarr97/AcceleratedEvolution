Code in Preprocessing was used to create the input files
Code in Preparing_For_Sampling was used to create the species and background probability matrices to assign substitutions to a lineage and gene during sampling.
Code in Sampling_And_Testing was used to do the sampling and then calculate p-values.
Code in Scripts_Making_Scripts is just helper python scripts to create .sh files and submit jobs.
Please see README and scripts in each folder for explanations of how things work.

For posthoc analysis, the jupyter notebooks were used as follows:
Cortical_development_analysis.ipynb was used for the analysis of HAG/CAG/unaccelerated gene enrichment for cortical marker genes.
Decreasing_CA.ipynb was used to do analysis and plotting for parts of figure 5 related to integration and validation with neuronal bulk ATAC-seq and ASE variance
MASH.ipynb was used to prepare the input files for mashr on ChromBPNet-related results, filter, and then make associated plots.  mash.R was used to actually run mashr.
Noncoding_ChromBPNet_CA.ipynb was used to do all other analysis related to ChromBPNet-derived results (parts of figure 4 and figure 5)
Noncoding_PhyloP.ipynb was used to do the analysis related to accelerated evolution for non-coding regions (figure 3)
Nonsynonymous.ipynb was used to do the analysis related to accelerated evolution for protein-coding regions (most of figure 1)
Plotting_Regions.ipynb was used to visualize accelerated evolution (parts of various figures)
Protein_Stability.ipynb was used to analyze the RaSP predictions and relationship to substitutions in conserved sites (part of figure 1)
Testing_LossOfConstraint.ipynb was used to analyze polymorphisms near substitutions in highly conserved sites (parts of figure 5)
UTR.ipynb was used to do the analysis related to accelerated evolution for UTRs (figure 2)
Western.ipynb was used to do the statistics on the western blots (parts of figure 1)

get_regions_to_plot.py was used to get the regions from the human gtf that we plotted
get_SuppTable9.py was used to pull data for Supplemental Table 9
check_constraint2.py was used to do the analysis of polymorphisms near substitutions in highly conserved sites
test_cons_enrichment.py was used to do the analysis of enrichment of substitutions in highly conserved sites in accessible regions in each cell type.