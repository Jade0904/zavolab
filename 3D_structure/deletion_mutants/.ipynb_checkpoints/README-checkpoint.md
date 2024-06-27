# About

This folder's analysis is based on the paper *Computational modeling and prediction of deletion mutants* (https://doi.org/10.1016/j.str.2023.04.005)[1]. Inspired by this paper, we are interested in testing the protein structure metrics predicted by two protocols: (1) using AlphaFold2 only; (2) using Alphafold2 + Rosetta Relax.


# Methods

In the paper, they use a small alpha-helical protein and biomolecular nuclear magnetic resonance (NMR) spectroscopy for generating a deletion mutant dataset that covers the entire amino acid sequence. To be specific, they use a small alpha-helical sterile alpha motif (SAM) domain. This original fasta file is stored under "fasta" folder. Here are the steps how I did my analysis.

First, use "generate_mutants.ipynb" to generate all the possible mutant sequences from the original sequence. This should result in 72 single deletion mutant sequences under the "fasta" folder.

Second, run AlphaFold2 on every mutant sequence. There is an example bash script at "bash_scripts/AlphaFold.sh". This should result in the AlphaFold prediction for every mutant under the specified folder. The process was run on clusters so it was submitted by batches using slurm.

Third, run Rosetta Relax on every AlphaFold result. The bash script using relax is at "bash_scripts/delRelax.sh". This should result in 3 outputs ending with ".out", ".sc" and ".log" under each AlphaFold folder. These are the Rosetta Relax outputs. I tried to use the same version of Rosetta as the original paper, which is not the latest version today, for consistency with their data.

Also, there are CSP (chemical shift perturbation) data from the authors of the paper, under the folder "del_csp_data". These data only applies for soluble mutants.


# Results

The first thing I checked is whether the Figure 3B and the lower right part of Figure 3D can be reproduced. This was done in "score.ipynb". The exact same results are hard to get but we can basically have the same conclusions from the original figures and the ones I reproduced. It is reasonable to state that, when pLDDT scores averaged across all residues are plotted with delta delta G measures from Rosetta, they show a clear separation between soluble and insoluble variants.[1] And it's generally true that AlphaFold2-RosettaRelax protocol captures the relationship between delta T_M and score insufficiency.

Further, we are interested in whether the two protocols we chose differs. To be specific, we would like to know whether it's useful to apply Rosetta Relax after AlphaFold2 when predicting deletion mutants. For each mutants, we calculate the following things:

(1) RMSD between wildtype and mutants, after AlphaFold2 only (af_rmsd);

(2) RMSD between wildtype and mutants, wildtype after AlphaFold2 only, mutants after AlphaFold2 + Rosetta Relax (relax_rmsd_wtAF). (not so meaningful because of the inconsistency)

(3) RMSD between wildtype and mutants, after AlphaFold2 + Rosetta Relax (relax_rmsd_wtRelax).

(4) delta G (Rosetta energy score).

(5) delta delta G and absolute delta delta G, basically the difference between delta G for mutants and delta G for wildtype.

# Reference

[1] Woods H, Schiano D L, Aguirre J I, et al. Computational modeling and prediction of deletion mutants[J]. Structure, 2023, 31(6): 713-723. e3.

[2] Pak M A, Markhieva K A, Novikova M S, et al. Using AlphaFold to predict the impact of single mutations on protein stability and function[J]. Plos one, 2023, 18(3): e0282689.