# About

This folder contains the analysis for the indel (insertion / deletions) examples selected from Protein Data Bank (PDB). This continues the analysis from "deletion mutants", but instead of single deletions, we tried to test the protein structure metrics predicted by two protocols: (1) using AlphaFold2 only; (2) using AlphaFold2 + Rosetta Relax, on longer indels (at least > 1 aa) within pairs of isoform and canonical sequences.

# Methods

Here I listed the steps that I used to select indel examples from Protein Data Bank.

(1) Performed Multiple Sequence Alignment (MSA) on 70% and 90% clusters on protein sequence data from Protein Data Bank. All the data can be downloaded based on the information on this website: [File Download Services](https://www.rcsb.org/docs/programmatic-access/file-download-services). The clusters used are under the folder "clusters". For the sake of efficiency, only clusters with more than 80 sequences and less than 500 sequences were performed MSA.

(2) Extracted the gaps with length longer than 8 residues, and not overlapping with the beginning or the ending 9 residues of alignments. This step as well as MSA was performed using the R script "msa.R".

(3) Manually selected the first-round candidates of indel isoforms and their corresponding canonical sequence in the form of pairs, and excluded the artifacts where gaps do not represent indel positions inside the sequences.

(4) Downloaded the structure files (PDB format) for the first-round candidates, and extracted the sequences from them. This step was done using the bash script "extract_clean_pdb.sh" under the folder "bash_scripts". This is because in the Protein Data Bank, the sequence data provided in fasta files may not be exactly the same as they appear in structural data in PDB format.

(5) Filtered out the artifacts and kept the pairs with real indels (> 1 aa). In the end 14 pairs of indel candidates (28 sequences in total) were selected. The records of selection process are in this sheet: [Google Sheets - insertion examples](https://docs.google.com/spreadsheets/d/1Rb4ThXejUMN_6HelHK5oNvu3q2fgHL9hM2RWapIfIxI/edit?usp=sharing). Please ask for access if needed.

(6) For each sequence, AlphaFold2 followed by Rosetta Relax was performed. The example scripts for AlphaFold2 and Rosetta Relax were "AF.sh" and "relax.sh" under the folder "bash_scripts". After this step, there are 6 structures for each pair, in which each sequence has 3 structures: the one originally downloaded from Protein Data Bank (marked as "original"), the one predicted by AlphaFold2 (marked as "af"), and the one predicted by Rosetta Relax based on AlphaFold2 prediction (marked as "relax").

(7) Canonical forms and isoforms were re-labelled based on pLDDT scores after performing AlphaFold2, and the structures with higher pLDDT scores were marked as "canonical" or "reference", while the ones with lower pLDDT scores were "isoform".

# Results

## Metrics

Similar to Deletion Mutants, we would like to calculate several metrics within each pair, and compare the two protocols. Most of the metrics were calculated per common residue. Common residues was defined using pairwise sequence alignment. The following metrics were calculated:

(1) RMSD between canonical form and isoform (original, af, relax).

(2) The difference of Rosetta energy score (delta delta G) between canonical form and isoform (original, af, relax).

(3) Rosetta energy score (delta G) for each structure (original, af, relax; for both canonical form and isoform).

(4) pLDDT per residue (af).

(5) Relative solvent accessibility (RSA) per residue (original, af, relax; for both canonical form and isoform).

All these metrics were saved as csv files, under the folder "pairs_csv".

## Correlation between different metrics

### Correlations between RMSD values

The correlation coefficients (kendall's tau) were calculated between the RMSD in AlphaFold2 structures (RMSD_af) and RMSD in original structures downloaded from Protein Data Bank (RMSD_original), and between RMSD in structures after AlphaFold2 + Rosetta Relax (RMSD_relax) and RMSD_original. They are marked as "corr_af" and "corr_relax", correspondingly, and saved in "rmsd.csv" under the folder "correlations_csv".

### Correlations between delta delta G values

The correlation coefficients (kendall's tau) were calculated between the delta delta G values in AlphaFold2 structures (ddG_af) and delta delta G in original structures downloaded from Protein Data Bank (ddG_original), and between delta delta G in structures after AlphaFold2 + Rosetta Relax (ddG_relax) and ddG_original). They are marked as "corr_af" and "corr_relax", correspondingly, and saved in "ddG.csv" under the folder "correlations_csv".

# Visualization

Some visualizations are in "plots_examples.ipynb". Some interesting conclusions might be:

(1) If combining all pairs together into a big dataframe (just as saved in "pairs_all.csv" under the folder "pairs_csv") and calculate the correlation coefficients between RMSD_af and RMSD_original, and between RMSD_relax and RMSD_original, we can see the former is 0.460, while the latter is 0.580, which indicates after using the protocol AlphaFold2 + Rosetta Relax, the RMSD values per residue can reflect the true distances to a greater extent compared to the protocol AlphaFold2.

(2) If printing RMSD corr_af and corr_relax one against another (from "rmsd.csv"), we can see that in most of the pairs (9 out of 14, the ones above the diagonal), corr_relax is greater than corr_af. The two groups were also plotted in a box plot, and we can observe a slight increase in corr_relax compared to corr_af.

(3) We calculated the absolute difference between RMSD_af and RMSD_original (marked as Diff_af), and between RMSD_relax and RMSD_original (marked as Diff_relax), which represent how RMSD after two protocols different from the true values. Per-residue pLDDT scores were binned by the distance of 10. These were plotted as violin plots and box plots. While the differences between Diff_relax and Diff_af are very close to 0 for residues with higher pLDDT scores, Diff_relax values are lower than Diff_af for residues with lower pLDDT scores, indicating that the use of Rosetta Relax can to some extend make up for the unconvinced regions of the structures predicted by AlphaFold2 and lead the structures closer to their true formations.

(4) If calculating the correlation coefficients between delta G and RSA, we can see that the overall coefficients are higher in AlphaFold2 structures than the true structures; and there is also an increasing after performing Rosetta Relax compared to AlphaFold2. If categorized all residues into different groups based on hydrophobicity, it is obvious that the correlation coefficients are higher in hydrophobic residues than hydrophilic ones.