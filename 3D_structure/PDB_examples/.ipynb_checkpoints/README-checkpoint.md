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

Similar to Deletion Mutants, we would like to calculate several metrics within each pair, and compare the two protocols. Most of the metrics were calculated per common residue. Common residues was defined using pairwise sequence alignment. The following metrics were calculated:

(1) RMSD between canonical form and isoform (original, af, relax).

(2) The difference of Rosetta energy score ($\Delta$$\Delta$G) between canonical form and isoform (original, af, relax).

(3) Rosetta energy score ($\Delta$G) for each structure (original, af, relax; for both canonical form and isoform.