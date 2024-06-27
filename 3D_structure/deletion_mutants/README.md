# About

This folder's analysis is based on the paper *Computational modeling and prediction of deletion mutants* (https://doi.org/10.1016/j.str.2023.04.005)[1]. Inspired by this paper, we are interested in testing the protein structure metrics predicted by two workflows: (1) using AlphaFold2 only; (2) using Alphafold2 + Rosetta Relax.


# Methods

In the paper, they use a small alpha-helical protein and biomolecular nuclear magnetic resonance (NMR) spectroscopy for generating a deletion mutant dataset that covers the entire amino acid sequence. To be specific, they use a small alpha-helical sterile alpha motif (SAM) domain. This original fasta file is stored under "fasta" folder. Here are the steps how I did my analysis.

First, use "generate_mutants.ipynb" to generate all the possible mutant sequences from the original sequence. This should result in 72 single deletion mutant sequences under the "fasta" folder.

Second, run AlphaFold2 on every mutant sequence. There is an example bash script at "bash_scripts/AlphaFold.sh". This should result in the AlphaFold prediction for every mutant under the specified folder. The process was run on clusters so it was submitted by batches using slurm.

Third, run Rosetta Relax on every AlphaFold result.


# Reference

[1] Woods H, Schiano D L, Aguirre J I, et al. Computational modeling and prediction of deletion mutants[J]. Structure, 2023, 31(6): 713-723. e3.

[2] Pak M A, Markhieva K A, Novikova M S, et al. Using AlphaFold to predict the impact of single mutations on protein stability and function[J]. Plos one, 2023, 18(3): e0282689.