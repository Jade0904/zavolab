# About

This folder is about predicting protein structures for IPA candidates and the corresponding reference sequences using the protocol AlphaFold2 + Rosetta Relax, and calculating some metrics out of them.

# Methods

The workflows are very similar to the other two folders. All IPA candidates were extracted from "0716_LIHC_novel_transcripts.csv" and saved as fasta files under the folder "aa_seq", using "generate_fasta.ipynb". It concludes 10 pairs of good IPA candidates and 5 pairs of bad IPA candidates, pair 12 and 14 were excluded from the later analysis because they have large difference between isoform and reference.

Basically, AlphaFold2 followed by Rosetta Relax was performed on each selected sequence. The example scripts to run AlphaFold2 as well as Rosetta Relax were "AF_1.sh" and "relax_1.sh" under "bash_scripts" folder. The structure with lowest overall delta G score was extracted using the bash script "extract_relax_pdbs.sh" under the folder "bash_scripts". For the extracted pdb files, delta G scores per residue were calculated using "score_relax_perResidue.sh" under the same folder.

The following metrics for each pair were calculated in "stability_metrics":

(1) RMSD values calculated between backbone carbons of common residues between reference and isoform.

(2) Overall delta G values for reference and isoform.

(3) 10% quantile of per-residue delta G values for reference and isoform, for both highest 10% and lowest 10%.

(4) Overall pLDDT values for reference and isoform.

(5) 10% quantile of per-residue pLDDT scores for reference and isoform, for both highest 10% and lowest 10%.

(6) Median RSA for reference and isoform.

(7) Median RSA of hydrophilic residues for reference and isoform.

(8) Median RSA of hydrophobic residues for reference and isoform.

(9) Median RSA of each residue type for reference and isoform.

The per-residue information for each pair was saved in "pairs_csv" folder. The above metrics were saved in "metrics.csv".