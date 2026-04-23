Now that we have our lineage and background matrices, we can just run simulate_new_humchp.py for UTR or non-coding
We match the parameters that were used for the matrix creation
The run_sim*.sh scripts are examples of what was used to run the simulate_humchp_new.py

For nonsyn, we use amino acids instead of trinucleotides so have different scripts.
simulate_nonsyn.py is the baseline, simulate_geneset_nonsyn.py is for gene sets, and simulate_geneset_nonsyn_rm3.py is for removing the top three genes per gene set from each species

compute_z_and_p.py is used to compute the p-value, total corrected p-value, and distribution p-values for the non-coding approach using the normal approximation.
